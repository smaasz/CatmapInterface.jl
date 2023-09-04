"""
    ratelaw_TS(prefactor, Gf_IS, Gf_TS, T, activprod)

Computes the rate of an elementary reaction that passes through a transition state.

The rate law is based on the Arrhenius relation.
The change in Gibbs free energy between the initial and transition state is the needed activation energy.
A `prefactor` and the product of the activities `activprod` complete the rate law.
"""
function ratelaw_TS(prefactor, Gf_IS, Gf_TS, T, activprod)
    @local_phconstants k_B
    prefactor * exp(-(Gf_TS - Gf_IS) / (k_B * T)) * activprod
end

"""
    compute_free_energies!(free_energies, catmap_params::CatmapParams, σ, ϕ_we, ϕ, local_pH) 

Compute the Gibbs free energies of all species specified in the `catmap_params` by applying the specified correction modes.
"""
function compute_free_energies!(free_energies, catmap_params::CatmapParams, σ, ϕ_we, ϕ, local_pH)
    (; gas_thermo_mode, adsorbate_thermo_mode, electrochemical_thermo_mode) = catmap_params

    for (s, (; formation_energy)) in catmap_params.species_list
        free_energies[s] += formation_energy
    end

    gas_thermo_correction!             = getfield(@__MODULE__, gas_thermo_mode)
    adsorbate_thermo_correction!       = getfield(@__MODULE__, adsorbate_thermo_mode)
    electrochemical_thermo_correction! = getfield(@__MODULE__, electrochemical_thermo_mode)

    gas_thermo_correction!(free_energies, catmap_params)
    adsorbate_thermo_correction!(free_energies, catmap_params)

    # electrochemical corrections
    electrochemical_thermo_correction!(free_energies, catmap_params, σ, ϕ_we, ϕ, local_pH)
    nothing
end


"""
    create_reaction_network(catmap_params::CatmapParams)

Create a [ReactionSystem](https://docs.sciml.ai/Catalyst/stable/api/catalyst_api/#Catalyst.ReactionSystem) for the specified microkinetic model.

For each elementary reaction the [ratelaw_TS](@ref) is used.
No separate rate equations for the solvent (= H₂O) and the active sites are created but their activities can be specified as parameters.
For ficitious gases (OH⁻ and H⁺) and adsorbates the activity coefficients are assumed to be 1.
The activity coefficients of the gaseous species can specified as parameters.
The thermodynamical corrections to the DFT-data of the formation energies are applied according to the specified modes.
New modes can be added by the user by adding a function with the same name to the module. 
"""
function create_reaction_network(catmap_params::CatmapParams)
    (; species_list, T) = catmap_params

    @parameters σ ϕ_we ϕ local_pH
    free_energies = Dict(zip(keys(species_list), fill(Num(0.0), length(species_list))))
    compute_free_energies!(free_energies, catmap_params::CatmapParams, σ, ϕ_we, ϕ, local_pH)

    @variables t
    vars        = Dict{String, Num}()
    activ_coefs = Dict{String, Num}()
    for (s, sp) in species_list
        if s =="H2O_g" || isa(sp, SiteSpecies) # site species and the solvent are not considered proper species
            as      = Symbol("a$s")
            vars[s] = first(@parameters $as)
        elseif (isa(sp, FictiousSpecies) && s ≠ "ele_g") || isa(sp, AdsorbateSpecies) # fictious species and adsorbates have no activity coeff
            ss          = Symbol(s)
            vars[s]     = first(@species $ss(t))
        elseif (isa(sp, GasSpecies) && s ≠  "H2O_g")
            ss              = Symbol(s)
            vars[s]         = first(@species $ss(t))
            gs              = Symbol("γ$s")
            activ_coefs[s]  = first(@parameters $gs)
        end
    end
    rxs = Reaction[]
    for ((; educts, products, tstate), prefactor) in zip(catmap_params.reactions, catmap_params.prefactors)
        Gf_IS = Num(0.0)
        Gf_FS = Num(0.0)
        Gf_TS = Num(0.0)
        es = Num[]
        αs = Int[]
        ps = Num[]
        βs = Int[]
        af = Num(1.0)
        ar = Num(1.0)
        for (educt, factor) in educts
            sp = species_list[educt]
            if (educt =="H2O_g" || isa(sp, SiteSpecies))
                af *= vars[educt]^factor
            elseif (isa(sp, FictiousSpecies) && educt ≠ "ele_g") # activity is assumed 1 b/c their influence is in rate constant
                push!(es, vars[educt])
                push!(αs, factor)
            elseif isa(sp, AdsorbateSpecies) # activity coefficients are assumed to be 1
                push!(es, vars[educt])
                push!(αs, factor)
                af *= vars[educt]^factor
            elseif (isa(sp, GasSpecies) && educt ≠ "H2O_g")
                push!(es, vars[educt])
                push!(αs, factor)
                af *= (activ_coefs[educt] * vars[educt])^factor
            end
            Gf_IS += factor * free_energies[educt]
        end
        for (product, factor) in products
            sp = species_list[product]
            if (product =="H2O_g" || isa(sp, SiteSpecies))
                ar *= vars[product]^factor
            elseif (isa(sp, FictiousSpecies) && product ≠ "ele_g") # activity is assumed 1 b/c their influence is in rate constant
                push!(ps, vars[product])
                push!(βs, factor)
            elseif isa(sp, AdsorbateSpecies) # activity coefficients are assumed to be 1
                push!(ps, vars[product])
                push!(βs, factor)
                ar *= vars[product]^factor
            elseif (isa(sp, GasSpecies) && product ≠ "H2O_g")
                push!(ps, vars[product])
                push!(βs, factor)
                ar *= (activ_coefs[product] * vars[product])^factor
            end
            Gf_FS += factor * free_energies[product]
        end
        Gf_TS = isnothing(tstate) ? max(Gf_IS, Gf_FS) : free_energies[tstate.symbol]
        rxn_f = Reaction(ratelaw_TS(prefactor, Gf_IS, Gf_TS, T, af), es, ps, αs, βs; only_use_rate=true)
        rxn_r = Reaction(ratelaw_TS(prefactor, Gf_FS, Gf_TS, T, ar), ps, es, βs, αs; only_use_rate=true)
        push!(rxs, rxn_f)
        push!(rxs, rxn_r)
    end
    ReactionSystem(rxs, t, name = :microkinetics)
end

"""
    generate_function(rn::ReactionSystem, dvs::Vector{Num}, ps::Vector{Num})

Generate a mutating function from a `ReactionSystem` that computes the concentration fluxes due to the reaction.
"""
function generate_function(rn::ReactionSystem, dvs::Vector{Num}, ps::Vector{Num})
    @assert Set(dvs) == Set(species(rn))
    @assert Set(ps)  == Set(parameters(rn))

    species_map = speciesmap(rn)
    sys = convert(ODESystem, rn; combinatoric_ratelaws=false)
    eqs = equations(sys)
    rhss = [-1 * eqs[species_map[dv]].rhs for dv in dvs] # multiply by -1 because the orientation assumed in VoronoiFVM physics functions

    u = map(x -> ModelingToolkit.time_varying_as_func(ModelingToolkit.value(x), sys), dvs)
    p = map(x -> ModelingToolkit.time_varying_as_func(ModelingToolkit.value(x), sys), ps)
    t = ModelingToolkit.get_iv(sys)

    pre, sol_states = ModelingToolkit.get_substitutions_and_solved_states(sys, no_postprocess = false)

    f_expr = build_function(rhss, u, p, t; postprocess_fbody = pre, states = sol_states)[2]
    drop_expr(@RuntimeGeneratedFunction(@__MODULE__, f_expr))
end