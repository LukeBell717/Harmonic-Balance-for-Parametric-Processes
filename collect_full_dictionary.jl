function collect_full_dictionary(main_folder, sweep_directory, class)

    # include scripts that will be run
    include(main_folder * "/scripts/create_all_dictionaries.jl")
    include(main_folder * "/scripts/reorganize_dictionaries.jl")
    
    # create and save all dictionaries specified by sweep_directory 
    flux_num = create_all_dictionaries(sweep_directory, class)

    # collect full dictionary with key values according to harmonic
    # coefficients ("u1", "v1", "u2", "v2", etc.); each key will have
    # a value that is a matrix of dimensions max_branch x ω_N x ϕ_p_N x Φ_N
    full_dict = reorganize_dictionaries(sweep_directory, flux_num)
    
    # save final_dict 
    sweep_directory_number_str = match(r"\d+$", sweep_directory).match
    sweep_directory_number = parse(Int, sweep_directory_number_str)
    HarmonicBalance.save(sweep_directory * "/full_dict" * string(sweep_directory_number), full_dict)
    
    return full_dict
end 