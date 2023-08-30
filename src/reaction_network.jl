function rateconstants(prefactor, Gf_IS, Gf_TS, T, actionprod)
    @local_phconstants k_B
    prefactor * exp(-(Gf_TS - Gf_IS) / (k_B * T)) * actionprod
end

function create_reaction_network(catmap_params::CatmapParams)
    (; species_list, T, potential_reference_scale, Upzc) = catmap_params

    gas_thermo_correction!             = getfield(@__MODULE__, catmap_params.gas_thermo_mode)
    adsorbate_thermo_correction!       = getfield(@__MODULE__, catmap_params.adsorbate_thermo_mode)
    electrochemical_thermo_correction! = getfield(@__MODULE__, catmap_params.electrochemical_thermo_mode)

    free_energies = Dict(zip(keys(species_list), fill(Num(0.0), length(species_list))))
    # thermo corrections
    gas_thermo_correction!(free_energies, species_list, T)
    adsorbate_thermo_correction!(free_energies, species_list, T)

    @parameters σ ϕ_we ϕ local_pH a[1:length(species_list)]
    # electrochemical corrections
    electrochemical_thermo_correction!(free_energies, species_list, σ, ϕ_we, ϕ, Upzc, local_pH, T; potential_reference_scale)

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
        rxn_f = Reaction(rateconstants(prefactor, Gf_IS, Gf_TS, T, af), es, ps, αs, βs; only_use_rate=true)
        rxn_r = Reaction(rateconstants(prefactor, Gf_FS, Gf_TS, T, ar), ps, es, βs, αs; only_use_rate=true)
        push!(rxs, rxn_f)
        push!(rxs, rxn_r)
    end
    ReactionSystem(rxs, t, name = :microkinetics)
end