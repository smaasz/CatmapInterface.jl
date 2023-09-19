
"""
$(TYPEDEF)

Abstract supertype of species.
"""
abstract type AbstractSpecies end


"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct FictiousSpecies <:AbstractSpecies
    """
    Name of the species
    """
    species_name::String
    """
    Formation energy of the species in joule per mole
    """
    formation_energy::Float64
    """
    Partial pressure of the species in pascal
    """
    pressure::Float64
    function FictiousSpecies(; species_name, formation_energy, pressure)
        if pressure < 0.0
            throw(DomainError("pressure must be nonnegative"))
        end
        new(species_name, formation_energy, pressure)
    end
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct GasSpecies <: AbstractSpecies
    """
    Name of the species
    """
    species_name::String
    """
    Formation energy of the species in joule per mole
    """
    formation_energy::Float64
    """
    Partial pressure of the species in pascal
    """
    pressure::Float64
    """
    Normal modes of vibration of the species in m⁻¹
    """
    frequencies::Vector{Float64}
    """
    Henry constant used to convert between gas and liquid phase
    """
    henry_const::Union{Float64, Missing}
    function GasSpecies(; species_name, formation_energy, pressure, frequencies, henry_const)
        if pressure < 0.0
            throw(DomainError("pressure must be nonnegative"))
        end
        if any(frequencies .<= 0.0)
            throw(DomainError("all frequencies must be positive"))
        end
        if !ismissing(henry_const) && henry_const <= 0.0
            throw(DomainError("Henry constant must be positive"))
        end
        new(species_name, formation_energy, pressure, frequencies, henry_const)
    end
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct AdsorbateSpecies <: AbstractSpecies
    """
    Name of the species
    """
    species_name::String
    """
    Formation energy of the species in joule per mole
    """
    formation_energy::Float64
    """
    Surface coverage of the species
    """
    coverage::Float64
    """
    Site the species adsorbs to
    """
    site::String
    """
    Name of the surface the species adsorbs to
    """
    surface_name::String
    """
    Normal modes of vibration of the species in m⁻¹
    """
    frequencies::Vector{Float64}
    """
    Parameters specifying the electrochemcial correction to the formation energy due to surface charges: ``G_f(U=0, σ) ≈ G_f(U=0,σ=0)+a\\cdot σ + b\\cdot σ^2``
    """
    sigma_params::@NamedTuple{a::Union{Nothing, Float64}, b::Union{Nothing, Float64}}
    """
    Parameter used in [`CatmapInterface.first_order_adsorbate_interaction`](@ref) to specify the formation energy's dependence on the coverage
    """
    self_interaction_param::Float64
    """
    Parameter used in [`CatmapInterface.first_order_adsorbate_interaction`](@ref) to specify the formation energy's dependence on the coverage of other adsorbates
    """
    cross_interaction_params::Dict{String, Float64}
    function AdsorbateSpecies(; species_name, formation_energy, coverage, site, surface_name, frequencies, sigma_params=(;a=nothing, b=nothing), self_interaction_param=0.0, cross_interaction_params=Dict{String, Float64}())
        if coverage < 0.0 || coverage > 1.0
            throw(DomainError("coverage must be between 0 and 1"))
        end
        if any(frequencies .<= 0.0)
            throw(DomainError("all frequencies must be positive"))
        end
        new(species_name, formation_energy, coverage, site, surface_name, frequencies, sigma_params, self_interaction_param, cross_interaction_params)
    end
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct TStateSpecies <: AbstractSpecies
    """
    Name of the species
    """
    species_name::String
    """
    Formation energy of the species in joule per mole
    """
    formation_energy::Float64
    """
    Surface coverage of the species
    """
    coverage::Float64
    """
    Site the species adsorbs to
    """
    site::String
    """
    Name of the surface the species adsorbs to
    """
    surface_name::String
    """
    Normal modes of vibration of the species in m⁻¹
    """
    frequencies::Vector{Float64}
    """
    Transfer coefficient of the reaction passing through the specified transition state
    """
    β::Float64
    """
    Specifying the list of (stable) species (with stoichiometric factor) the transition state is between
    """
    between_species::Vector{Pair{String, Int}}
    """
    Parameters specifying the electrochemcial correction to the formation energy due to surface charges: ``G_f(U=0, σ) ≈ G_f(U=0,σ=0)+a\\cdot σ + b\\cdot σ^2``
    """
    sigma_params::@NamedTuple{a::Union{Nothing, Float64}, b::Union{Nothing, Float64}}
    """
    Parameter used in [`CatmapInterface.first_order_adsorbate_interaction`](@ref) to specify the formation energy's dependence on the coverage of other adsorbates
    """
    cross_interaction_params::Dict{String, Float64}
    function TStateSpecies(; species_name, formation_energy, coverage, site, surface_name, frequencies, β, between_species, sigma_params=(;a=nothing, b=nothing), cross_interaction_params=Dict{String, Float64}())
        if coverage < 0.0 || coverage > 1.0
            throw(DomainError("coverage must be between 0 and 1"))
        end
        if any(frequencies .<= 0.0)
            throw(DomainError("all frequencies must be positive"))
        end
        new(species_name, formation_energy, coverage, site, surface_name, frequencies, β, between_species, sigma_params, cross_interaction_params)
    end
end
#TStateSpecies(; formation_energy, coverage, site, surface_name, frequencies, sigma_params::Vector{Float64}) = TStateSpecies(; formation_energy, coverage, site, surface_name, frequencies, sigma_params=(; a=sigma_params[1], b=sigma_params[2]))

@kwdef struct InteractionResponseParams
    slope::Float64      = 1.0
    cutoff::Float64     = 0.25
    smoothing::Float64  = 0.05
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
@kwdef struct SiteSpecies <: AbstractSpecies
    """
    Formation energy of the site in joule per mole
    """
    formation_energy::Float64
    """
    Name of the site
    """
    site_name::String
    """
    Parameters used to specify the exact form of the `interaction_response_function` in [`CatmapInterface.first_order_adsorbate_interaction`](@ref)
    """
    interaction_response_params::InteractionResponseParams = InteractionResponseParams()
end