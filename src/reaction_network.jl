function rateconstants(prefactor, Gf_IS, Gf_TS, T, activity_coeffs)
    @local_phconstants k_B
    prefactor * exp(-(Gf_TS - Gf_IS) / (k_B * T)) * prod(activity_coeffs)
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
        ss  = Symbol(s)
        if (isa(sp, GasSpecies) && s ≠ "H2O_g") || isa(sp, AdsorbateSpecies) || (isa(sp, FictiousSpecies) && s ≠ "ele_g")
            vars[s] = first(@species $ss(t))
        end
        if isa(sp, GasSpecies) || s == "ele_g" || isa(sp, SiteSpecies)
            as              = Symbol("γ$s")
            activ_coefs[s]  = first(@parameters $as)
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
        af = Num[]
        ar = Num[]
        for (educt, factor) in educts
            sp = species_list[educt]
            if (isa(sp, GasSpecies) && educt ≠ "H2O_g") || isa(sp, AdsorbateSpecies) || (isa(sp, FictiousSpecies) && educt ≠ "ele_g")
                push!(es, vars[educt])
                push!(αs, factor)
            end
            if isa(sp, GasSpecies) || educt == "ele_g" || isa(sp, SiteSpecies)
                push!(af, activ_coefs[educt]^factor)
            end
            Gf_IS += factor * free_energies[educt]
        end
        for (product, factor) in products
            sp = species_list[product]
            if (isa(sp, GasSpecies) && product ≠ "H2O_g") || isa(sp, AdsorbateSpecies) || (isa(sp, FictiousSpecies) && product ≠ "ele_g")
                push!(ps, vars[product])
                push!(βs, factor)
            end
            if isa(sp, GasSpecies) || product == "ele_g" || isa(sp, SiteSpecies)
                push!(ar, activ_coefs[product]^factor)
            end
            Gf_FS += factor * free_energies[product]
        end
        Gf_TS = isnothing(tstate) ? max(Gf_IS, Gf_FS) : free_energies[tstate.symbol]
        rxn_f = Reaction(rateconstants(prefactor, Gf_IS, Gf_TS, T, af), es, ps, αs, βs)
        rxn_r = Reaction(rateconstants(prefactor, Gf_FS, Gf_TS, T, ar), ps, es, βs, αs)
        push!(rxs, rxn_f)
        push!(rxs, rxn_r)
    end
    ReactionSystem(rxs, t, name = :microkinetics)
end