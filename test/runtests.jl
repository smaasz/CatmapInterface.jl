using CatmapInterface
using Test

catmap_params   = parse_catmap_input("catmap_CO2R_template.mkm")
rn              = create_reaction_network(catmap_params)

@testset "CatmapInterface.jl" begin
    @test length(reactions(rn)) = 5
    @test length(species(rn)) == 11
end

# θ       = [0.2, 0.8]
# prob    = SteadyStateProblem(kinetic_model!, θ, p = (rxn_parameter_grid[77,1], rxn_parameter_grid[77,2], [1.0, 1.0, 1.0])) 
# sol     = solve(prob, DynamicSS(Rodas5()))

# println(descriptor_grid[77, :])
# println(sol)

#@testset "CatmapInterface.jl" begin
    # Write your tests here.
#end
