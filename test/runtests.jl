using CatmapInterface
using Test
using DifferentialEquations
using CSV

function get_catmap_coverage_map(csvfile_path)

    file = CSV.File(csvfile_path, header = ["d1", "d2", "c1", "c2"])
    entries = []
    for row in file
        push!(entries, ([row.d1, row.d2], [row.c1, row.c2]))
    end

    Dict(entries)
end


(descriptor_grid, rxn_parameter_grid, kinetic_model)   = get_kinetic_model("./CO_oxidation.mkm", "./energies.txt")

kinetic_model! = convert_to_julia(kinetic_model)

catmap_coverage_map = get_catmap_coverage_map("./coverage_map_catmap.csv")
coverage_map        = compute_coverage_map(descriptor_grid, rxn_parameter_grid, kinetic_model!)

@testset "CatmapInterface.jl" for descriptor_values in eachrow(descriptor_grid)
    @test catmap_coverage_map[descriptor_values] == coverage_map[descriptor_values]
end

# θ       = [0.2, 0.8]
# prob    = SteadyStateProblem(kinetic_model!, θ, p = (rxn_parameter_grid[77,1], rxn_parameter_grid[77,2], [1.0, 1.0, 1.0])) 
# sol     = solve(prob, DynamicSS(Rodas5()))

# println(descriptor_grid[77, :])
# println(sol)

#@testset "CatmapInterface.jl" begin
    # Write your tests here.
#end
