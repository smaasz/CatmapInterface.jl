using .CatmapInterface
using Test
using DifferentialEquations

(grid, kinetic_model)   = get_kinetic_model("/home/steffen/Documents/masterarbeit/catmap-1/tutorials/2-creating_microkinetic_model/CO_oxidation.mkm", "/home/steffen/Documents/masterarbeit/catmap-1/tutorials/2-creating_microkinetic_model/energies.txt")
kinetic_model!          = convert_to_julia(kinetic_model)


θ       = [0.4, 0.4]
prob    = SteadyStateProblem(kinetic_model!, θ, p = (grid[1,1], grid[1,2], [1.0, 1.0, 1.0])) 
sol     = solve(prob, DynamicSS(QNDF()))
println(sol)

#@testset "CatmapInterface.jl" begin
    # Write your tests here.
#end
