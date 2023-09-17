using CatmapInterface
using Test
using Catalyst
using DelimitedFiles
using PyCall

if haskey(ENV, "CATMAP")
    pushfirst!(pyimport("sys")."path", ENV["CATMAP"])
end

const eV = 1.602176634e-19

# The following python code block defines a function to extract the free energies of all reactants given a CatMAP input file
py"""
from catmap import ReactionModel
import catmap

def catmap_kinetic_model(setup_file, theta):

    # ReactionModel is the main class that is initialized with a setup-file
    model = ReactionModel(setup_file=setup_file)

    # some solver parameters have to be set manually (?!)
    import mpmath as mp
    model.solver._mpfloat = mp.mpf
    model.solver._math = mp
    model.solver._matrix = mp.matrix

    model.solver.compile() # compiles all templates, here (rate_constants) are needed

    # Set up interaction model.
    if model.adsorbate_interaction_model == 'first_order':
        interaction_model = \
            catmap.thermodynamics.FirstOrderInteractions(model)
        interaction_model.get_interaction_info()
        response_func = interaction_model.interaction_response_function
        if not callable(response_func):
            int_function = getattr(interaction_model,
                                    response_func+'_response')
            interaction_model.interaction_response_function = int_function
        model.thermodynamics.__dict__['adsorbate_interactions'] = interaction_model
    elif model.adsorbate_interaction_model in ['ideal',None]:
        model.thermodynamics.adsorbate_interactions = None
    else:
        raise AttributeError(
                'Invalid adsorbate_interaction_model specified.')

    descriptor_values = [descriptor_range[0] for descriptor_range in model.descriptor_ranges]
    rxn_parameters = model.scaler.get_rxn_parameters(descriptor_values)
    n_tot = len(model.adsorbate_names) + len(model.transition_state_names)
    energies = rxn_parameters[:n_tot]
    if len(rxn_parameters) == n_tot + n_tot**2:
        interaction_vector = rxn_parameters[-n_tot**2:]
    elif len(rxn_parameters) == n_tot:
        interaction_vector = [0]*n_tot**2
    F = model.interaction_response_function if hasattr(model, 'interaction_response_function') else None
    theta = [theta[species] for species in model.adsorbate_names] + len(model.transition_state_names) * [0.0]
    _, Gf, _ = model.interaction_function(theta, energies, interaction_vector, F)

    energies = {
        k:v for k, v in zip(
            model.adsorbate_names + model.transition_state_names,
            Gf
        )
    }
    for k, v in zip(model.gas_names, model.solver._gas_energies):
        energies[k] = v
    for k, v in zip(model.site_names, model.solver._site_energies):
        energies[k] = v

    return energies
"""

"""
    instantiate_catmap_template(template_file_path, params)

Instantiate a template file by inserting the parameters in the `NamedTuple` `params`.
"""
function instantiate_catmap_template(template_file_path, params)
    (; σ, ϕ_we, ϕ, local_pH, ϕ_pzc, T) = params
    input_instance_string = open(template_file_path, "r") do template_file
        read(template_file, String)
    end

	replacements = Dict(
		r"descriptor_ranges.?=.*" =>"descriptor_ranges = [[$ϕ_we, $ϕ_we], [$T, $T]]",
		r"voltage_diff_drop.?=.*" => "voltage_diff_drop = $ϕ",
		r"pH.?=.*" => "pH = $local_pH",
		r"Upzc.?=.*" => "Upzc = $ϕ_pzc",
		r"\nsigma_input.?=.*" => "\nsigma_input = $σ/0.01", # in μF/cm^2
	)
    input_instance_string = replace(input_instance_string, replacements...)

    input_instance_file_name = joinpath(dirname(template_file_path), "CO2R.mkm")
    open(input_instance_file_name, "w") do input_instance_file
        write(input_instance_file, input_instance_string)
    end

    return input_instance_file_name
end

"""
    compute_catmap_free_energies(catmap_instance_path, params)

Compute the free energies of all reactants in the microkinetic model from CatMAP.
"""
function compute_catmap_free_energies(catmap_instance_path, params)
    (; θ) = params
    starting_dir = pwd()
    catmap_instance_path = abspath(catmap_instance_path)
    cd(dirname(catmap_instance_path))
    try
        free_energies = py"catmap_kinetic_model"(catmap_instance_path, θ)
    finally
        cd(starting_dir)
    end
 	free_energies = convert(Dict{String, Float64}, free_energies)
    delete!(free_energies, "g")
    for (k, v) in free_energies
        if length(k) == 1
            delete!(free_energies, k)
            free_energies["_$k"] = v
        end
    end
    free_energies
end

"""
    compute_interface_free_energies(catmap_instance_path, params)

Compute the free energies of all reactants in the microkinetic model from the CatmapInterface.
"""
function compute_interface_free_energies(catmap_instance_path, params)
    (; θ, σ, ϕ_we, ϕ, local_pH) = params
    catmap_params   = parse_catmap_input(catmap_instance_path)    
    free_energies   = Dict(zip(keys(catmap_params.species_list), fill(0.0, length(catmap_params.species_list))))
    CatmapInterface.compute_free_energies!(free_energies, catmap_params, θ, σ, ϕ_we, ϕ, local_pH)
    free_energies
end

"""
    test_free_energies(catmap_template_path, params; rtol=1.0e-5)

Create a testset where the free energies of all reactants in the microkinetic model computed by CatMAP and the CatmapInterface are compared.
"""
function test_free_energies(catmap_template_path, params; rtol=1.0e-5)
    catmap_instance_path    = instantiate_catmap_template(catmap_template_path, params)
    catmap_free_energies    = compute_catmap_free_energies(catmap_instance_path, params)
    interface_free_energies = compute_interface_free_energies(catmap_instance_path, params)
    rm(catmap_instance_path)
    @assert keys(catmap_free_energies) == keys(interface_free_energies) "The names of the species don't match: $(keys(catmap_free_energies)) vs. $(keys(interface_free_energies))"
    @testset "species=$species" for species in keys(interface_free_energies)
        @test isapprox(catmap_free_energies[species], interface_free_energies[species]/eV; rtol)
    end 
end


const Cgap = 0.2 # in F/m^2
models = [
    (;  
        model                   = "CO₂-Reduction on Au with hbond corrections", 
        catmap_template_path    = "../data/Au-model-hbond/catmap_CO2R_template.mkm", 
        params_set              = map(
            row -> (; zip([:θ                                                                                                                                                               ,:ϕ_we      ,:local_pH  ,:T     ,:ϕ_pzc     ,:ϕ     ,:σ     ], [row; Cgap * (row[2] - row[6] - row[5])])...),
            eachrow([
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.80
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.72
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.64
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.56
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.48
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.40
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.32
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.16
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.08
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.24
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        +0.00
            ])
        )
    ),
    (;  
        model                   = "CO₂-Reduction on Au with simple corrections", 
        catmap_template_path    = "../data/Au-model-simple/catmap_CO2R_template.mkm", 
        params_set              = map(
            row -> (; zip([:θ                                                                                                                                                               ,:ϕ_we      ,:local_pH  ,:T     ,:ϕ_pzc     ,:ϕ     ,:σ     ], [row; Cgap * (row[2] - row[6] - row[5])])...),
            eachrow([
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.80
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.72
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.64
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.56
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.48
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.40
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.32
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.16
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.08
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.24
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        +0.00
            ])
        )
    ),
    (;  
        model                   = "CO₂-Reduction on Cu with simple corrections", 
        catmap_template_path    = "../data/Liu-model-simple/catmap_CO2R_template.mkm", 
        params_set              = map(
            row -> (; zip([:θ                                                                                                                                                               ,:ϕ_we      ,:local_pH  ,:T     ,:ϕ_pzc     ,:ϕ     ,:σ     ], [row; Cgap * (row[2] - row[6] - row[5])])...),
            eachrow([
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.80
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.72
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.64
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.56
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.48
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.40
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.32
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.16
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.08
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.24
                            Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        +0.00
            ])
        )
    ),
    (;  
        model                   = "CO₂-Reduction on Cu with simple corrections and first order adsorbate corrections", 
        catmap_template_path    = "../data/Liu-model-first-order/catmap_CO2R_template.mkm", 
        params_set              = map(
            row -> (; zip([:θ                                                                                                                                                               ,:ϕ_we      ,:local_pH  ,:T     ,:ϕ_pzc     ,:ϕ     ,:σ     ], [row; Cgap * (row[2] - row[6] - row[5])])...),
            eachrow([
                            Dict("OCCO_t" => 0.2, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.80
                            Dict("OCCO_t" => 0.2, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.72
                            Dict("OCCO_t" => 0.2, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.64
                            Dict("OCCO_t" => 0.2, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.56
                            Dict("OCCO_t" => 0.2, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.48
                            Dict("OCCO_t" => 0.2, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.40
                            Dict("OCCO_t" => 0.2, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.32
                            Dict("OCCO_t" => 0.2, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.16
                            Dict("OCCO_t" => 0.2, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.08
                            Dict("OCCO_t" => 0.2, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.24
                            Dict("OCCO_t" => 0.2, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        +0.00
            ])
        )
    ),
]

@testset "CatmapInterface.jl" begin
    @testset "$model" for (model, catmap_template_path, params_set) in models
        @testset "$(params.ϕ_we)" for params in params_set
            test_free_energies(catmap_template_path, params)
        end
    end
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
