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