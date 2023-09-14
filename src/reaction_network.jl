"""
$(SIGNATURES)

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
$(SIGNATURES)

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
$(SIGNATURES)

Create a [ReactionSystem](https://docs.sciml.ai/Catalyst/stable/api/catalyst_api/#Catalyst.ReactionSystem) for the specified microkinetic model.

For each elementary reaction the [`CatmapInterface.ratelaw_TS`](@ref) is used.
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

    function process_reaction_side(reactants)
        @local_unitfactors mol dm
        Gf = Num(0.0)
        rs = Num[]
        γs = Int[]
        a  = Num(1.0)
        for (reactant, factor) in reactants
            sp = species_list[reactant]
            if (reactant =="H2O_g" || isa(sp, SiteSpecies))
                a *= vars[reactant]^factor
            elseif (isa(sp, FictiousSpecies) && reactant ≠ "ele_g") # activity is assumed 1 b/c their influence is in rate constant
                push!(rs, vars[reactant])
                push!(γs, factor)
            elseif isa(sp, AdsorbateSpecies) # activity coefficients are assumed to be 1
                push!(rs, vars[reactant])
                push!(γs, factor)
                a *= (vars[reactant])^factor
            elseif (isa(sp, GasSpecies) && reactant ≠ "H2O_g")
                push!(rs, vars[reactant])
                push!(γs, factor)
                a *= (activ_coefs[reactant] * vars[reactant])^factor
            end
            Gf += factor * free_energies[reactant]
        end
        Gf, rs, γs, a
    end

    rxs = Reaction[]
    for ((; educts, products, tstate), prefactor) in zip(catmap_params.reactions, catmap_params.prefactors)
        (Gf_IS, es, αs, af) = process_reaction_side(educts)
        (Gf_FS, ps, βs, ar) = process_reaction_side(products)
        Gf_TS = isnothing(tstate) ? max(Gf_IS, Gf_FS) : mapreduce(x->free_energies[first(x)]^last(x), +, tstate.components) #free_energies[tstate.name]
        rxn_f = Reaction(ratelaw_TS(prefactor, Gf_IS, Gf_TS, T, af), es, ps, αs, βs; only_use_rate=true)
        rxn_r = Reaction(ratelaw_TS(prefactor, Gf_FS, Gf_TS, T, ar), ps, es, βs, αs; only_use_rate=true)
        push!(rxs, rxn_f)
        push!(rxs, rxn_r)
    end
    ReactionSystem(rxs, t, name = :microkinetics)
end

"""
$(SIGNATURES)

Transform the (micro-)kinetic model from surface/gas-reactions to surface/electrolyte-reactions using Henry's law.
"""
function liquidize(odesys::ODESystem, catmap_params::CatmapParams)
    @local_unitfactors bar
    (; species_list) = catmap_params

    sts     = states(odesys)
    ps      = parameters(odesys)

    usubs = Pair{SymbolicUtils.BasicSymbolic{Real}, SymbolicUtils.BasicSymbolic{Real}}[]
    csubs = Pair{Num, Num}[]
    psubs = Pair{SymbolicUtils.BasicSymbolic{Real}, SymbolicUtils.BasicSymbolic{Real}}[]
    @variables t
    for st in sts
        sp = species_list[string(Symbolics.operation(Symbolics.value(st)))]
        if isa(sp, GasSpecies)
            (; henry_const) = sp
            ss          = Symbol("$(sp.species_name)_aq")
            var         = first(@variables $ss(t))
            push!(usubs, st => Symbolics.value(var))
            push!(csubs, st => var / henry_const / bar)

            gs          = Symbol("γ$(sp.species_name)_aq")
            activ_coef  = first(@parameters $gs)
            p = getproperty(odesys, Symbol("γ$(sp.species_name)_g"); namespace=false)
            p = Symbolics.value(p)
            push!(psubs, p => Symbolics.value(activ_coef))
        end
    end

    new_eqs = Equation[]
    for eq in equations(odesys)
        lhs = expand_derivatives(substitute(eq.lhs, Dict(usubs)))
        rhs = substitute(eq.rhs, Dict(csubs..., psubs...))
        push!(new_eqs, Equation(lhs, rhs))
    end
    structural_simplify(ODESystem(new_eqs, t, replace(sts, usubs...), replace(ps, psubs...); name=odesys.name))
end

"""
$(SIGNATURES)

Generate a mutating function from a `ReactionSystem` that computes the concentration fluxes due to the reaction.
"""
function generate_function(rn::ReactionSystem; dvs::Vector{Tval}=species(rn), ps::Vector{Tval}=parameters(rn)) where {Tval <: Union{SymbolicUtils.BasicSymbolic{Real}, Num}}
    @assert Set(dvs) == Set(species(rn))
    @assert Set(ps)  == Set(parameters(rn))

    species_map = speciesmap(rn)
    sys = convert(ODESystem, rn; combinatoric_ratelaws=false)
    eqs = equations(sys)
    rhss = [-1 * eqs[species_map[dv]].rhs for dv in dvs] # multiply by -1 because the orientation assumed in VoronoiFVM physics functions

    u = map(x -> ModelingToolkit.time_varying_as_func(Symbolics.value(x), sys), dvs)
    p = map(x -> ModelingToolkit.time_varying_as_func(Symbolics.value(x), sys), ps)
    t = ModelingToolkit.get_iv(sys)

    pre, sol_states = ModelingToolkit.get_substitutions_and_solved_states(sys, no_postprocess = false)

    f_expr = build_function(rhss, u, p, t; postprocess_fbody = pre, states = sol_states)[2]
    drop_expr(@RuntimeGeneratedFunction(@__MODULE__, f_expr))
end

"""
$(SIGNATURES)

Generate a mutating function from a `ODESystem` that computes the concentration fluxes due to the reaction.
"""
function generate_function(sys::ODESystem; dvs=states(sys), ps=parameters(sys))
    @assert Set(dvs) == Set(states(sys))
    @assert Set(ps)  == Set(parameters(sys))


    #state_map = Dict(zip(states(sys), length(states(sys))))
    state_map = Dict([st => i for (i, st) in enumerate(states(sys))])
    eqs = equations(sys)
    rhss = [-1 * eqs[state_map[dv]].rhs for dv in dvs] # multiply by -1 because the orientation assumed in VoronoiFVM physics functions

    u = map(x -> ModelingToolkit.time_varying_as_func(Symbolics.value(x), sys), dvs)
    p = map(x -> ModelingToolkit.time_varying_as_func(Symbolics.value(x), sys), ps)
    t = ModelingToolkit.get_iv(sys)

    pre, sol_states = ModelingToolkit.get_substitutions_and_solved_states(sys, no_postprocess = false)

    f_expr = build_function(rhss, u, p, t; postprocess_fbody = pre, states = sol_states)[2]
    drop_expr(@RuntimeGeneratedFunction(@__MODULE__, f_expr))
end