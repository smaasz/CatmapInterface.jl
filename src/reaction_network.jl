function rateconstants(prefactor, Gf_IS, Gf_TS, T, activity_coeffs)
    @local_phconstants k_B
    prefactor * exp(-(Gf_TS - Gf_IS) / (k_B * T)) * prod(activity_coeffs)
end

function create_reaction_network(catmap_params::CatmapParams)
    (; species_list, T, potential_reference_scale, Upzc) = catmap_params

    species_strings = keys(species_list)
    species_strings_to_idxs = Dict(zip(species_strings, 1:length(species_strings)))

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

    species_symbols = Expr[]
    for s in species_strings
        ss = Symbol(s)
        push!(species_symbols, :($ss(t)))
    end
    @variables t
    species_num = eval(:(@species $(species_symbols...)))
    rxs         = Reaction[]
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
            push!(es, species_num[species_strings_to_idxs[educt]])
            push!(αs, factor)
            push!(af, a[species_strings_to_idxs[educt]])
            Gf_IS += free_energies[educt]
        end
        for (product, factor) in products
            push!(ps, species_num[species_strings_to_idxs[product]])
            push!(βs, factor)
            push!(ar, a[species_strings_to_idxs[product]])
            Gf_FS += free_energies[product]
        end
        Gf_TS = isnothing(tstate) ? max(Gf_IS, Gf_FS) : free_energies[tstate.symbol]
        rxn_f = Reaction(rateconstants(prefactor, Gf_IS, Gf_TS, T, af), es, ps, αs, βs)
        rxn_r = Reaction(rateconstants(prefactor, Gf_FS, Gf_TS, T, ar), ps, es, βs, αs)
        push!(rxs, rxn_f)
        push!(rxs, rxn_r)
    end
    ReactionSystem(rxs, t, name = :microkinetics)
end