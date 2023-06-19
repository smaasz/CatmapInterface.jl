function compute_coverage_map(descriptor_grid, rxn_parameter_grid, kinetic_model!, θ_init, p; use_odesolver = false)

    entries = []

    for (descriptor_values, rxn_parameters) in ProgressBar(zip(eachrow(descriptor_grid), eachrow(rxn_parameter_grid)))
    
        prob    = SteadyStateProblem(kinetic_model!, θ_init, p = (rxn_parameters[1], rxn_parameters[2], p))
        
        if use_odesolver
            sol     = solve(prob, DynamicSS(Rodas5()))
        else #use NonlinearSolve.jl
            #prob    = NonlinearProblem(prob)
            #sol     = solve(prob, NewtonRaphson())
            sol     = solve(prob, SSRootfind())
        end

        # add to coverage map
        push!(entries, (descriptor_values, sol))

        θ_init  = sol

    end

    Dict(entries)
end