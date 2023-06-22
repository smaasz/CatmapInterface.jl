function convert_to_julia(steady_state_expressions) 
    steady_state_expressions = replace(steady_state_expressions, r"(\w+?)\[(\d*?)\]" => s"\1[\2 + 1]")
    steady_state_expressions = replace(steady_state_expressions, r"theta" => s"θ")
    steady_state_expressions = replace(steady_state_expressions,  r"mpf\('(.*)'\)" => s"\1") # r"mpf\('(.*)'\)" => s"BigFloat('\1')"
    steady_state_expressions = replace(steady_state_expressions,  r"(\h)\[(.*?)\]\h*\*\h*(.*?)([\h|\n])" => s"\1\2 * ones(T, \3)\4")
    # r"([^a-zA-Z])len\((.*?)\)" => s"\1length(\2)"
    return steady_state_expressions 
end

function convert_pythonfloat(x)
    return Float64(BigFloat(pycall(py"str", String, x)))
end

function get_coverage_equations(steady_state_expressions)
    steady_state_system_template = """
    function kinetic_model!(dθ_dt, θ::AbstractVector{T}, params, t) where T
        
        kf, kr, p   = params 
        r           = zeros(T, length(kf))
        
        $steady_state_expressions
    end
    """
    return eval(Meta.parse(steady_state_system_template))
end

function get_elementary_rates_equations(steady_state_expressions)
    elementary_rates_template = """
    function elementary_rates(kf, kr, θ::AbstractVector{T}, p) where T
        r       = zeros(T, length(kf))
        dθ_dt   = zeros(T, length(θ))
        
        $steady_state_expressions
        
        return r    
    end
    """
    return eval(Meta.parse(elementary_rates_template))
end

function get_catmap_output(setup_file_path, energies_file_path)
    
    setup_file_path         = realpath(setup_file_path)
    energies_file_path      = realpath(energies_file_path)
    (~, setup_file_name)    = splitdir(setup_file_path)
    (~, energies_file_name) = splitdir(energies_file_path)

    workingdir = pwd()
    tmpdir = mktempdir()
    cd(tmpdir)
    tmp_setup_file_path     = cp(setup_file_path, joinpath(tmpdir, setup_file_name))
    tmp_energies_file_path  = cp(energies_file_path, joinpath(tmpdir, energies_file_name))

    descriptor_grid, rxn_parameter_grid, steady_state_expressions, elementary_rxns, gas_names = py"catmap_kinetic_model"(tmp_setup_file_path)

    rxn_parameter_grid  = ((x) -> convert_pythonfloat.(x)).(rxn_parameter_grid)
    descriptor_grid     = ((x) -> convert_pythonfloat.(x)).(descriptor_grid)


    steady_state_expressions    = convert_to_julia(steady_state_expressions)
    kinetic_model!              = get_coverage_equations(steady_state_expressions)
    elementary_rates            = get_elementary_rates_equations(steady_state_expressions)

    function compute_rates(rxn_parameters, θ, p)
       
        kf = rxn_parameters[1]
        kr = rxn_parameters[2]
    
        rs =  elementary_rates(kf, kr, θ, p)

        return rs
    end

    cd(workingdir)

    return descriptor_grid, rxn_parameter_grid, kinetic_model!, compute_rates
end