function get_kinetic_model(setup_file_path, energies_file_path)
    
    setup_file_path         = realpath(setup_file_path)
    energies_file_path      = realpath(energies_file_path)
    (~, setup_file_name)    = splitdir(setup_file_path)
    (~, energies_file_name) = splitdir(energies_file_path)

    workingdir = pwd()
    tmpdir = mktempdir()
    cd(tmpdir)
    tmp_setup_file_path     = cp(setup_file_path, joinpath(tmpdir, setup_file_name))
    tmp_energies_file_path  = cp(energies_file_path, joinpath(tmpdir, energies_file_name))

    grid, kin_model = py"catmap_kinetic_model"(tmp_setup_file_path)

    grid = ((x) -> convert_pythonfloat.(x)).(grid)

    cd(workingdir)

    return grid, kin_model

end

function convert_pythonfloat(x)
    return Float64(BigFloat(pycall(py"str", String, x)))
end

function convert_to_julia(kin_model)

    # verify that the kinetic model fulfills the assumptions

    # remove last three lines
    kin_model = join(split(kin_model, "\n")[1:end-5], "\n") * "end"
    
    #kin_model = replace(kin_model, r"def (?<functionname>.+?)\((?<arguments>.*)\):" => s"function \g<functionname>!(dtheta_dt, \g<arguments>)")
    kin_model = replace(kin_model, r".*def.*" => s"function kinetic_model!(dtheta_dt::AbstractVector{T}, theta::AbstractVector{T}, ps, t) where T\n\n    kf, kr, p = ps\n")
    
    kin_model = replace(kin_model, r"dtheta_dt = .*" => s"")
    
    kin_model = replace(kin_model, r"(\w+?)\[(\d*?)\]" => s"\1[\2 + 1]")
    kin_model = replace(kin_model, r"mpf\('(.*)'\)" => s"\1")
    #kin_model = replace(kin_model, r"mpf\('(.*)'\)" => s"BigFloat('\1')")
    kin_model = replace(kin_model, r"(\h)\[(.*?)\]\h*\*\h*(.*?)([\h|\n])" => s"\1\2 * ones(T, \3)\4")
    kin_model = replace(kin_model, r"([^a-zA-Z])len\((.*?)\)" => s"\1length(\2)")

    #return kin_model
    return eval(Meta.parse(kin_model))
end