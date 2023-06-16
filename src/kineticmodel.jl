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
    # do things
    println("test")
    println(tmp_setup_file_path)
    a = py"catmap_kinetic_model"(tmp_setup_file_path)

    cd(workingdir)

    return a

end

function convert_to_julia(kin_model)
    
    kin_model = replace(kin_model, r"def (.+):" => s"function \1") * "end"
    kin_model = replace(kin_model, r"(\w+?)\[(\d*?)\]" => s"\1[\2 + 1]")
    kin_model = replace(kin_model, r"mpf\('(.*)'\)" => s"\1")
    kin_model = replace(kin_model, r"(\h)\[(.*?)\]\h*\*\h*(.*?)([\h|\n])" => s"\1\2 * ones(\3)\4")
    kin_model = replace(kin_model, r"([^a-zA-Z])len\((.*?)\)" => s"\1length(\2)")

    return kin_model
end