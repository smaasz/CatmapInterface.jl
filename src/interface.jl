struct Species
    energy::Float64
    frequencies::Vector{Float64}
    pressure::Float64
    sigma_params::@NamedTuple{a::Float64, b::Float64}
    function Species(; pressure, sigma_params)
        if pressure < 0.0
            throw(DomainError("pressure must be nonnegative"))
        end
        if any(frequencies .<= 0.0)
            throw(DomainError("all frequencies must be positive"))
        end
        new(energy, frequencies, pressure, sigma_params)
    end
end

"""
Main data structure for interfacing CatMAP
"""
struct CatmapParams
    rn::ReactionSystem
    prefactors::Vector{Float64}
    species_definitions::Dict{Num, Species}
    gas_thermo_mode::Symbol
    adsorbate_thermo_mode::Symbol
    electrochemical_thermo_mode::Symbol
    function CatmapParams(; rn, prefactors, species, gas_thermo_mode, adsorbate_thermo_mode, electrochemical_thermo_mode)       
        if !(length(prefactors == numreactions(rn)))
            throw(ArgumentError("The number of prefactors must match the number of reactions"))
        end
        if any(prefactors .<= 0.0)
            throw(DomainError("all prefectors must be positive"))
        end
        # check that for each species in the reaction network parameters are provided
        ks = keys(species)
        function isspecified(sp)
            any(map(isequal(sp), keys(species)))
        end
        for sp in species(rn)
            if !isspecified(sp) 
                throw(ArgumentError("$sp not specified in the argument species_definitions"))
            end
        end
        # check that the provided modes are defined
        for mode in [gas_thermo_mode, adsorbate_thermo_mode, electrochemical_thermo_mode]
            if !isdefined(@__MODULE__, mode)
                throw(ArgumentError("$(String(:mode))=$mode is not implemented"))
            end
        end
        new(rn, prefactors, species_definitions, gas_thermo_mode, adsorbate_thermo_mode, electrochemical_thermo_mode)
    end
end