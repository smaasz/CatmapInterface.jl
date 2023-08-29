using Test
using CatmapInterface
using Catalyst

catmap_params   = parse_catmap_input("catmap_CO2R_template.mkm")
rn              = create_reaction_network(catmap_params)

@testset "CatmapInterface.jl" begin
    @test numreactions(rn) == 2*5
    @test numspecies(rn) == 10
end

# θ       = [0.2, 0.8]
# prob    = SteadyStateProblem(kinetic_model!, θ, p = (rxn_parameter_grid[77,1], rxn_parameter_grid[77,2], [1.0, 1.0, 1.0])) 
# sol     = solve(prob, DynamicSS(Rodas5()))

# println(descriptor_grid[77, :])
# println(sol)

#@testset "CatmapInterface.jl" begin
    # Write your tests here.
#end
