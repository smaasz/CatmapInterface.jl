"""
"""
function ideal_adsorbate_interaction(energies, catmap_params, coverages)
    nothing
end

    
# c_0 = 0.0
# if θ_tot <= x0
#     c_0 = 0.0
# elseif θ_tot <= x1
#     α = slope/(2*(x1-x0))
#     c_0 = (α * (θ_tot - x0)^2)/θ_tot
#     #dC = alpha*(1-(x0/theta_tot)**2)
#     #d2C = (2*alpha*x0**2)/(theta_tot**3)
# else
#     c_0 = slope*(θ_tot - cutoff)/θ_tot
#     # dC = slope*(cutoff/(theta_tot**2))
#     # d2C = (-2*slope*cutoff)/(theta_tot**3)
# end
# c_0

function smooth_piecewise_linear_interaction_response(θ_tot, interaction_response_params::InteractionResponseParams)
    (; slope, cutoff, smoothing) = interaction_response_params
    θ_1 = min(max(0.0, θ_tot - (cutoff - smoothing)), 2 * smoothing)
    θ_2 = max(0.0, θ_tot - (cutoff + smoothing))

    α   = smoothing > 0 ? slope / (4 * smoothing) : 0.0
    (α * θ_1^2 + slope * θ_2)
end

function linear_interaction_response(θ_tot, interaction_response_params::InteractionResponseParams)
    (; slope) = interaction_response_params
    smooth_piecewise_linear_interaction_response(θ_tot, InteractionResponseParams(; slope, cutoff=0.0, smoothing=0.0))
end

function piecewise_linear_interaction_response(θ_tot, interaction_response_params::InteractionResponseParams)
    (; slope, cutoff) = interaction_response_params
    smooth_piecewise_linear_interaction_response(θ_tot, InteractionResponseParams(; slope, cutoff, smoothing=0.0))
end

function geometric_mean_cross_interaction(ϵ_s, ϵ_os)
    √(ϵ_s * ϵ_os)
end

function arithemtic_mean_cross_interaction(ϵ_s, ϵ_os)
    (ϵ_s + ϵ_os) / 2
end

function neglect_cross_interaction(ϵ_s, ϵ_os)
    0.0
end

function intermediate_state_transition_cross_interaction(ϵ_rs)
    ϵ_rs * 0.5
end

function neglect_transition_cross_interaction(ϵ_rs)
    ϵ_rs * 0.5
end

"""
$(SIGNATURES)

Extract the interaction term if needed by applying the `cross_interaction_function`
"""
function get_interaction_term(s::String, sp::AbstractSpecies, os::String, osp::AbstractSpecies, cross_interaction_function)
    ϵ = nothing
    if s == os
        ϵ = sp.self_interaction_param
    else
        ϵ = get(sp.cross_interaction_params, os, nothing)
        if isnothing(ϵ)
            ϵ = get(osp.cross_interaction_params, s, nothing)
        end
        if isnothing(ϵ)
            ϵ_s  = sp.self_interaction_param
            ϵ_os = osp.self_interaction_param
            ϵ    = cross_interaction_function(ϵ_s, ϵ_os)
        end
    end 
    return ϵ
end

function coverage_of_site(site, species_list, θ)
    θ_tot = 0.0
    for (s, sp) in species_list
        if isa(sp, AdsorbateSpecies) && (site == sp.site)
            θ_tot += θ[s]
        end
    end
    θ_tot
end

"""
$(SIGNATURES)

Add correction terms to the adsorbation energies based on first order interactions between the adsorbates. 

This adsorbation interaction model is expained in ![CatMAP's documentation](https://catmap.readthedocs.io/en/latest/topics/including_adsorbate_adsorbate_interactions.html#coverage-dependent-adsorption-eneriges) and in this ![issue](https://github.com/smaasz/CatmapInterface.jl/issues/10)
"""
function first_order_adsorbate_interaction(energies, catmap_params::CatmapParams, θ)
    @local_unitfactors eV
    (; species_list, adsorbate_interaction_params) = catmap_params
    (; interaction_response_function, cross_interaction_mode, transition_state_cross_interaction_mode) = adsorbate_interaction_params

    response_function                           = getfield(@__MODULE__, Symbol(interaction_response_function, "_interaction_response"))
    cross_interaction_function                  = getfield(@__MODULE__, Symbol(cross_interaction_mode, "_cross_interaction"))
    transition_state_cross_interaction_function = getfield(@__MODULE__, Symbol(transition_state_cross_interaction_mode, "_transition_cross_interaction"))

    for (s, sp) in species_list
        if isa(sp, AdsorbateSpecies)
            (; site)                        = sp
            θ_tot                           = coverage_of_site(site, species_list, θ)
            (; interaction_response_params) = species_list["_$site"]
            response_value                  = response_function(θ_tot, interaction_response_params)
            for (os, osp) in species_list
                if isa(osp, AdsorbateSpecies)
                    ϵ = get_interaction_term(s, sp, os, osp, cross_interaction_function)
                    sp.cross_interaction_params[os] = ϵ
                    osp.cross_interaction_params[s] = ϵ
                    energies[s] += response_value * ϵ * ((θ[os] + 1.0e-15)/(θ_tot + 1.0e-15)) * eV
                end
            end
        end
    end

    for (s, sp) in species_list
        if isa(sp, TStateSpecies)
            (; site)                        = sp
            θ_tot                           = coverage_of_site(site, species_list, θ)
            (; interaction_response_params) = species_list["_$site"]
            response_value                  = response_function(θ_tot, interaction_response_params)
            for (os, osp) in species_list
                if isa(osp, AdsorbateSpecies)
                    ϵ = get(sp.cross_interaction_params, os, nothing)
                    if isnothing(ϵ)
                        ϵ = get(osp.cross_interaction_params, s, nothing)
                    end
                    if isnothing(ϵ)
                        ϵ_rs = mapreduce(+, sp.between_species) do reactant
                            if isa(species_list[reactant], AdsorbateSpecies)
                                return species_list[reactant].cross_interaction_params[os]
                            else
                                return 0.0
                            end
                        end
                        ϵ = transition_state_cross_interaction_function(ϵ_rs)
                    end
                    sp.cross_interaction_params[os] = ϵ
                    osp.cross_interaction_params[s] = ϵ
                    energies[s] += response_value * ϵ * ((θ[os] + 1.0e-15)/(θ_tot + 1.0e-15)) * eV
                end 
            end
        end
    end
end

"""
$(SIGNATURES)

Add thermodynamic correction terms for all gas species using the ideal gas approximation.
"""
function ideal_gas(energies::Dict{String, Tval}, catmap_params::CatmapParams) where Tval <: Real
    @local_unitfactors eV
    @local_phconstants h c_0
    
    ideal_gas_params = py"ideal_gas_params"
    (; species_list, T) = catmap_params
    for (s, sp) in species_list
        if isa(sp, GasSpecies)
            (; species_name, frequencies) = sp
            try
                (symmetrynumber, geometry, spin) = ideal_gas_params[s]
                energies[s] += py"get_thermal_correction_ideal_gas"(T, frequencies * (h * c_0 / eV), symmetrynumber, geometry, spin, species_name) * eV
            catch e
                if isa(e, KeyError)
                    throw(ArgumentError("$s has no specified ideal gas parameters"))
                else
                    throw(e)
                end
            end
        end
    end
end

"""
$(SIGNATURES) 

Add thermodynamic correction terms for all adsorbed species using the harmonic adsorbate approximation.
"""
function harmonic_adsorbate(energies, catmap_params::CatmapParams)
    @local_unitfactors eV
    @local_phconstants h c_0
    (; species_list, T) = catmap_params
    for (s, sp) in species_list
        if isa(sp, AdsorbateSpecies)
            (; frequencies) = sp
            energies[s] += py"get_thermal_correction_adsorbate"(T, frequencies * (h * c_0 / eV)) * eV
        end
    end
    for (s, sp) in species_list
        
        if isa(sp, TStateSpecies)
            if !isempty(sp.frequencies)
                (; frequencies) = sp
                energies[s] += py"get_thermal_correction_adsorbate"(T, frequencies * (h * c_0 / eV)) * eV
            else
                (; between_species) = sp
                for bs in between_species
                    if isa(species_list[bs], AdsorbateSpecies)
                        thermo_correction = energies[bs]
                        energies[s] += 0.5 * thermo_correction
                    end
                end
            end
        end

    end
    nothing
end

"""
$(SIGNATURES) 

Add electrochemical pH-correction terms for all influenced species.
"""
function _get_echem_corrections(energies, catmap_params::CatmapParams, σ, ϕ_we, ϕ, local_pH)
    @local_unitfactors eV
    (; species_list, potential_reference_scale, T) = catmap_params
    if haskey(species_list, "H_g") || haskey(species_list, "OH_g")
        G_H2 = energies["H2_g"]
        if potential_reference_scale == "SHE"
            G_H = G_H2 * 0.5 - (0.0592 * local_pH / 298.14 * T) * eV
        elseif potential_reference_scale == "RHE"
            G_H = G_H2 * 0.5
        end
        G_H2O = energies["H2O_g"]
        G_OH = G_H2O - G_H
        if haskey(species_list, "H_g")
            energies["H_g"] += G_H
        end
        if haskey(species_list, "OH_g")
            energies["OH_g"] += G_OH
        end
    end
    nothing
end

"""
$(SIGNATURES) 

Add electrochemical correction terms for all influenced species using ...
"""
function simple_electrochemical(energies, catmap_params::CatmapParams, σ, ϕ_we, ϕ, local_pH)
    @local_unitfactors eV
    (; species_list, Uref) = catmap_params
    # simple_electrochem_corrections
    if haskey(energies, "ele_g")
        energies["ele_g"] += -(ϕ_we - ϕ) * eV
    end
    for (s, sp) in species_list
        if isa(sp, TStateSpecies) && occursin("ele", sp.species_name)
            (; β) = sp
            energies[s] += (-(ϕ_we - ϕ) + β * (ϕ_we - ϕ - Uref)) * eV
        end
    end
    # pH_correction
    _get_echem_corrections(energies, catmap_params, σ, ϕ_we, ϕ, local_pH)
    nothing
end

"""
$(SIGNATURES) 

Add electrochemical correction terms for all influenced species using ...
"""
function hbond_surface_charge_density(energies, catmap_params::CatmapParams, σ, ϕ_we, ϕ, local_pH)
    @local_unitfactors eV
    # simple_electrochem
    simple_electrochemical(energies, catmap_params, σ , ϕ_we, ϕ, local_pH)
    hbond_dict = py"hbond_dict"
    (; species_list) = catmap_params
    # hbond_electrochemical
    for (s, sp) in species_list
        if isa(sp, AdsorbateSpecies)
            (; species_name) = sp
            if haskey(hbond_dict, species_name)
                energies[s] += hbond_dict[species_name] * eV
            end
        end
    end
    # hbond_surface_charge_density
    for (s, sp) in species_list
        if (isa(sp, AdsorbateSpecies) || isa(sp, TStateSpecies))
            (; a, b) = sp.sigma_params
            energies[s] += (a * σ + b * σ^2) * eV
        end
    end
    nothing
end