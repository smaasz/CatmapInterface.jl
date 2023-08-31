using CatmapInterface
using Test
using Catalyst
using DelimitedFiles

const eV = 1.602176634e-19

function instantiate_catmap_template(template_file_path, σ, ϕ_we, ϕ, local_pH; ϕ_pzc=0.16, T=298.0)
    input_instance_string = open(template_file_path, "r") do template_file
        read(template_file, String)
    end

	replacements = Dict(
		r"descriptor_ranges.?=.*" =>"descriptor_ranges = [[$ϕ_we, $ϕ_we], [$T, $T]]",
		r"voltage_diff_drop.?=.*" => "voltage_diff_drop = $ϕ",
		r"pH.?=.*" => "pH = $local_pH",
		r"Upzc.?=.*" => "Upzc = $ϕ_pzc",
		r"\nsigma_input.?=.*" => "\nsigma_input = $σ",
	)
    input_instance_string = replace(input_instance_string, replacements...)

    input_instance_file_name = "CO2R.mkm"
    open(input_instance_file_name, "w") do input_instance_file
        write(input_instance_file, input_instance_string)
    end

    return input_instance_file_name
end

function compare_catmap_output(catmap_output_path; rtol=1.0e-5)
    (dc, dh) = readdlm(catmap_output_path, ',', Float64; header=true)
    
    function findparams(row)
        ϕ_we        = row[findfirst(isequal("ϕ_we"), dh[1,:])]
        ϕ           = row[findfirst(isequal("ϕ"), dh[1,:])]
        ϕ_pzc       = row[findfirst(isequal("ϕ_pzc"), dh[1,:])]
        local_pH    = row[findfirst(isequal("local_pH"), dh[1,:])]
        T           = row[findfirst(isequal("T"), dh[1,:])]
        σ           = row[findfirst(isequal("σ"), dh[1,:])]
        (; ϕ_we, ϕ, ϕ_pzc, local_pH, T, σ)
    end

    function test_catmap_output(row)
        (; ϕ_we, ϕ, ϕ_pzc, local_pH, T, σ) = findparams(row)

        input_instance_file = instantiate_catmap_template("catmap_CO2R_template.mkm", σ, ϕ_we, ϕ, local_pH; ϕ_pzc, T)
        catmap_params       = parse_catmap_input(input_instance_file)
        rm(input_instance_file)
        
        free_energies = Dict(zip(keys(catmap_params.species_list), fill(0.0, length(catmap_params.species_list))))
        CatmapInterface.compute_free_energies!(free_energies, catmap_params, σ, ϕ_we, ϕ, local_pH)

        @testset "species=$s" for s in collect(keys(free_energies))
            catmap_output = 0.0
            ss = startswith(s, "_") ? s[2:end] : s
            catmap_output = row[findfirst(isequal(ss), dh[1,:])]
            @test isapprox(catmap_output, free_energies[s] / eV; rtol)
        end
    end

    @testset "$(findparams(row))" for row in eachrow(dc)
        test_catmap_output(row)
    end
end

catmap_params   = parse_catmap_input("catmap_CO2R_template.mkm")
rn              = create_reaction_network(catmap_params)

@testset "CatmapInterface.jl" begin
    @test numreactions(rn) == 2*4
    @test numspecies(rn) == 6
    compare_catmap_output("./catmap_data.csv")
end



# md"""
# ## CatMAP
# """

# # ╔═╡ 1f020e47-e4cb-4202-80f7-39421e4ce471
# function instantiate_catmap_template(template_file_path, ϕ_we, ϕ, local_pH)
    
#     input_instance_string = open(template_file_path, "r") do template_file
#         read(template_file, String)
#     end

# 	replacements = Dict(
# 		r"descriptor_ranges.?=.*" =>"descriptor_ranges = [[$ϕ_we, $ϕ_we], [$T, $T]]",
# 		r"voltage_diff_drop.?=.*" => "voltage_diff_drop = $ϕ",
# 		r"pH.?=.*" => "pH = $local_pH",
# 		r"Upzc.?=.*" => "Upzc = $ϕ_pzc",
# 		r"\nsigma_input.?=.*" => "\nsigma_input = $((ϕ_we - ϕ - ϕ_pzc) * C_gap/(μF/cm^2))",
# 	)
#     input_instance_string = replace(input_instance_string, replacements...)

#     input_instance_file_name = "CO2R.mkm"
#     open(input_instance_file_name, "w") do input_instance_file
#         write(input_instance_file, input_instance_string)
#     end

#     return input_instance_file_name
# end

# # ╔═╡ 9d4a2505-9c24-4ed2-bf02-f6598f296f11
# begin
# 	template_file_path = abspath("../catmap_CO2R_template.mkm")
# 	isfile(template_file_path)
# end

# # ╔═╡ 6120ec04-80c8-423c-a92a-44583ba9bac7
# catmap_df = let
# py"""
# from catmap import ReactionModel

# def catmap_kinetic_model(setup_file):

# 	# ReactionModel is the main class that is initialized with a setup-file
# 	model = ReactionModel(setup_file=setup_file)

# 	# some solver parameters have to be set manually (?!)
# 	import mpmath as mp
# 	model.solver._mpfloat = mp.mpf
# 	model.solver._math = mp
# 	model.solver._matrix = mp.matrix

# 	model.solver.compile() # compiles all templates, here (rate_constants) are needed

# 	# define descriptor grid
# 	## model.resolution can either be a list of integers or an integer, if it's only
# 	## an integer, then duplicate the same resolution for all descriptors to normalize
# 	## model.resolution to a list
# 	if isinstance(model.resolution, int):
# 		model.resolution = [model.resolution] * len(model.descriptor_names)

# 	ndescriptor_tuples = 1
# 	for res in model.resolution:
# 		ndescriptor_tuples *= res

# 	## Note: define linspace to keep dependencies minimal
# 	def linspace(start, stop, nitems):
# 		if start >= stop:
# 			return [start]
# 		delta = (stop - start)/(nitems-1)
# 		r = []
# 		for i in range(nitems):
# 			r.append(start + i * delta)
# 		return r

# 	## list of the ranges of the descriptors
# 	descriptor_ranges = [linspace(_range[0], _range[1], res) for _range, res in zip(model.descriptor_ranges, model.resolution)]

# 	## The combinations of the descriptors are ordered in lexicographical order
# 	## This function returns the position of the combination at position idx in the
# 	## descriptor grid
# 	def get_idx_tuple(idx, lengths):
# 		idx_tuple = []
# 		rest = idx
# 		for length in lengths:
# 			idx_tuple.append(rest % length)
# 			rest = rest // length
# 		return idx_tuple
	
# 	descriptor_grid = [0] * ndescriptor_tuples
# 	rxn_params_grid = [0] * ndescriptor_tuples
	
# 	for i in range(ndescriptor_tuples):		
# 		descriptor_values = [
# 			descriptor_range[idx] for idx, descriptor_range in zip(get_idx_tuple(i, model.resolution), descriptor_ranges)
# 		]
# 		descriptor_grid[i] = descriptor_values
# 		rxn_params_grid[i] = {
# 			k:v for k, v in zip(
# 				model.adsorbate_names + model.transition_state_names,
# 				model.scaler.get_rxn_parameters(descriptor_values)
# 			)
# 		}
# 		for k, v in zip(model.gas_names, model.solver._gas_energies):
# 			rxn_params_grid[i][k] = v
# 		for k, v in zip(model.site_names, model.solver._site_energies):
# 			rxn_params_grid[i][k] = v
# 	return rxn_params_grid[0], model.species_definitions
# """
# 	df = nothing
# 	for ϕ in range(ϕ_we, 0.0, 11)
# 		σ = (ϕ_we - ϕ - ϕ_pzc) * C_gap
# 		params = Dict("ϕ_we" => ϕ_we, "ϕ" => ϕ, "σ" => σ)
		
# 		catmap_instance_path = instantiate_catmap_template(template_file_path, ϕ_we, ϕ, local_pH)
# 		catmap_instance_path = abspath(joinpath(".", catmap_instance_path))
# 		energies, specd = py"catmap_kinetic_model"(catmap_instance_path)
# 		energies = convert(Dict{String, Float64}, energies)
		
# 		if isnothing(df)
# 			df = DataFrame(params..., energies...)
# 		else
# 			push!(df, merge(params, energies))
# 		end
# 	end
# 	df
# end




# θ       = [0.2, 0.8]
# prob    = SteadyStateProblem(kinetic_model!, θ, p = (rxn_parameter_grid[77,1], rxn_parameter_grid[77,2], [1.0, 1.0, 1.0])) 
# sol     = solve(prob, DynamicSS(Rodas5()))

# println(descriptor_grid[77, :])
# println(sol)

#@testset "CatmapInterface.jl" begin
    # Write your tests here.
#end
