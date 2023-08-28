function create_reaction_network(catmap_params::CatmapParams)
    species_list = catmap_params.species_list

    gas_thermo_correction!             = getfield(@__MODULE__, catmap_params.gas_thermo_mode)
    adsorbate_thermo_correction!       = getfield(@__MODULE__, catmap_params.adsorbate_thermo_mode)
    electrochemical_thermo_correction! = getfield(@__MODULE__, catmap_params.electrochemical_thermo_mode)

    free_energies = Dict(zip(keys(species_list), fill(Num(0.0), length(species_list))))
    # thermo corrections
    gas_thermo_correction!(free_energies, species_list)
    adsorbate_thermo_correction!(free_energies, species_list)

    @parameters σ ϕ_we
    @variables t
    @species 
    # electrochemical corrections
    electrochemical_thermo_correction!(free_energies, species_list, σ, ϕ_we)

    end
end