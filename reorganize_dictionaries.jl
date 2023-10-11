# ==============================================================================
# function assemble all dictionaries from a given flux sweep into one dictionary
# ==============================================================================



function assemble_dictionaries(sweep_directory, flux_num)
    
    # create empty dictionary whose keys are the harmonic coefficients of the ansatz
    data = HarmonicBalance.load(sweep_directory * "/results1.jld2")
    variable_list = [string(element) for element in data.problem.variables]
    flux_dict = Dict{String, Vector{Any}}(variable => [] for variable in variable_list)
    
    # loop through all flux points 
    for dict_num in 1:(flux_num + 1)
        
        try
            # load dictionary for each flux point from directory specified by "sweep_directory" 
            dictionary = HarmonicBalance.load(sweep_directory * "/dict$dict_num.jld2")
            
            # loop through data for each key in dict and append data to flux_dict; 
            # this creates another "index" for flux_dict, which serve as the index 
            # for a given flux point
            for key in keys(dictionary)
                variable = dictionary[key]
                push!(flux_dict[key], variable)
            end
        catch
            continue
        end
    end
    return flux_dict
end



# =======================================================================================
# function to find the maximum number of branches for any given flux point
# =======================================================================================



function find_max_branches(flux_dict, flux_num)
    
    # create a list that will store the max number of branches for each harmonic coefficient
    all_branch_list = []
    
    # loop through all harmonic coefficients in flux_dict
    for key in keys(flux_dict)
        # empty list to store numbers of all branches for a given flux point
        branch_list = []
            # loop through all flux points 
            for dict_num in 1:flux_num
                # find the number of branches for a given flux point and append this to branch_list
                num_branch = size(flux_dict[key][dict_num])[1]
                push!(branch_list, num_branch)
            end 
        # find max branch number for a given harmonic coefficient and appned this to all_branch_list
        max_branch = maximum(branch_list)
        push!(all_branch_list, max_branch)
    end 
    
    # check whether the max number of branches for each harmonic coefficient is the same
    if all(x -> x == all_branch_list[1], all_branch_list)
        branch_max = all_branch_list[1]
    else 
        println("There is a mismatch in the number of branches")
    end
    
    return branch_max
    
end



# ========================================================================================
# function to make the dictionary associated to each flux point the same number of entries
# ========================================================================================



function homogenize_dictionaries(flux_dict, flux_num, branch_max)
    
    # create an empty dictionary with the harmonic coefficients as keys
    flux_dict_homogenized = Dict{String, Any}(key => [] for key in keys(flux_dict))
        
    # loop through all flux points 
    for dict_num in 1:flux_num

        # for each flux point, loop through all harmonic coefficients in the ansatz
        for key in keys(flux_dict)

            # get variable from dictionary; flux_dict has dimensions flux_point x branch_number x frequency x power; 
            # flux_dict[key][dict_num] has diensions branch_number x frequency x power
            variable = flux_dict[key][dict_num]

            # get the number of branches, frequency points and power points for a given flux point 
            branch_num = size(variable)[1]
            frequency_num = size(variable)[2] 
            power_num = size(variable)[3]
                
            # append matrix of NaNs to original variable so that all variables 
            # will have the same number of "branches"; 
            # NaN matrix has dimensions append_branch x frequency_num x power_num
            branch_append = branch_max - branch_num 
            NaN_matrix = fill(NaN, branch_append, frequency_num, power_num); 
            filled_matrix = vcat(variable, NaN_matrix); 
            # append to flux_dict_4d
            push!(flux_dict_homogenized[key], filled_matrix)
        end 
    end
    
    # collect number of points in sweep of frequency and power
    f_N = size(flux_dict_homogenized["u1"][1])[2]
    E_p_N = size(flux_dict_homogenized["u1"][1])[3]
    
    # return reshaped dictionary
    return f_N, E_p_N, flux_dict_homogenized
    
end



# =======================================================================================================================
# function to reshape the dictionary so that the final dictionary has dimensions flux x branch_number x frequency x power
# =======================================================================================================================



function reshape_dictionaries(flux_dict_homogenized, branch_max, f_N, E_p_N, flux_num)
    
    # loop through all keys in flux_dict_homogenized
    for key in keys(flux_dict_homogenized)
        
        # flatten and reshape dictionary so that it has dimensions 
        # branch_number x frequency x power x flux_points
        flat_flux_dict = collect(Iterators.flatten(flux_dict_homogenized[key]))
        reshaped_flux_dict = reshape(flat_flux_dict, (branch_max, f_N, E_p_N, flux_num))
        flux_dict_homogenized[key] = reshaped_flux_dict
    end 
    
    return flux_dict_homogenized
end



# ==============================================================================
# function assemble all dictionaries from a given flux sweep into one dictionary
# ==============================================================================



function reorganize_dictionaries(sweep_directory, flux_num)
    
    # assemble all dictionaries from a flux sweep into one dictionary
    flux_dict = assemble_dictionaries(sweep_directory, flux_num)
    
    # find the maximum number of branches for any given flux point
    branch_max = find_max_branches(flux_dict, flux_num)
    
    # prepare all dictionaries to be reshaped, and find the number of frequency points and power points 
    f_N, E_p_N, flux_dict_homogenized = homogenize_dictionaries(flux_dict, flux_num, branch_max)
    
    # reshape dictionareis so that each key in the dictionary has dimensions 
    # branch_max x frequency_num x power_num x flux_num
    final_dict = reshape_dictionaries(flux_dict_homogenized, branch_max, f_N, E_p_N, flux_num)
    
    return final_dict
end