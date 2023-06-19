function compute_coverage_map(descriptor_grid, rxn_parameter_grid, kinetic_model!; p = [1.0, 1.0, 1.0])

    entries = []

    θ_init  = [0.5, 0.5]
    for (descriptor_values, rxn_parameters) in zip(descriptor_grid, rxn_parameter_grid)
    
        prob    = SteadyStateProblem(kinetic_model!, θ_init, p = (rxn_parameters[1], rxn_parameters[2], p))
        sol     = solve(prob, DynamicSS(Rodas5()))

        # add to coverage map
        push!(entries, (descriptor_values, sol))

        θ_init  = sol

    end

    Dict(entries)
end