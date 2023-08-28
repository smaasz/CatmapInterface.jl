function ideal_gas(energies, species_list)
    for k in keys(energies)
        energies[k] += 1.0
    end
end

function harmonic_adsorbate(energies, species_list)
    for k in keys(energies)
        energies[k] += 1.0
    end
end

function hbond_surface_charge_density(energies, species_list, σ, ϕ_we, local_pH, T)
    for k in keys(energies)
        energies[k] += σ
    end
end