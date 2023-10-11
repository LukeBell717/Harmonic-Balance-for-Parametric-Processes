using HarmonicBalance



# ==========================================
# function to filter data by applying a mask 
# ==========================================



function filter_data(data, variable, class, not_class)
    
    # create array from 1 to the total number of branches
    all_branches = 1:HarmonicBalance.branch_count(data) 
    
    # collect variable across all branches
    variable_raw = transform_solutions(data, variable, branches = all_branches)
    
    # apply mask
    variable_masked = HarmonicBalance._apply_mask(variable_raw, HarmonicBalance._get_mask(
        data, class, not_class, branches = all_branches))
    
    
    # collect number of branches and size of sweep
    num_branches = length(variable_masked[1, 1])
    sweep_size = size(variable_masked)
    
    # reshape data 
    shape = (num_branches, sweep_size...) # specify size of matrix into which we will reshape data
    variable_masked_flatten = collect(Iterators.flatten(variable_masked)) # flatten data
    variable_reshaped = real(reshape(variable_masked_flatten, shape)) # reshape data
    
        
    # get mask to delete all rows that contain only NaNs
        # all(isnan.(matrix), dims=2) evaluates  whether the entries of the rows of the matrix are equal to NaN
        # .! elementwise negates every boolean returned by all(isnan.(matrix), dims=2)
        # since .!(all(isnan.(matrix), dims=2)) is a 6 x 1 BitMatrix, [:, 1] converts this to a BitVector
    row_mask = (.!(all(isnan.(variable_reshaped), dims=(2,3))))[:, 1, 1]

    # get new filtered matrix; this is the data for one harmonic coefficient
    variable_filtered = variable_reshaped[row_mask, :, :]
    
    return variable_filtered

end



# ============================================================================
# function to assemble dictionary for all harmonic variables in a given ansatz
# ============================================================================



function create_dictionary(sweep_directory, result_num, class)
    
    # load data; note this loads data from one simulation in a given sweep
    data = HarmonicBalance.load(sweep_directory * "/results"*string(result_num)*".jld2")
    
    # classiy vacuum solutions for signal and intermediate photons
    classify_solutions!(data, "sqrt(u1^2 + v1^2) < 1e-2", "pump_vacuum")
    classify_solutions!(data, "sqrt(u2^2 + v2^2) < 1e-2", "signal_vacuum")
    classify_solutions!(data, "sqrt(u3^2 + v3^2) < 1e-2", "intermediate_vacuum")
    
    # collect a list of all the harmonic coefficients in the ansatz (u1, v1, ... , un, vn)
    variable_list = [string(element) for element in data.problem.variables]; 

    # create a dictionary where each entry in variable_list is the key for data processed by filter_data
    dictionary = Dict() 
    for variable in variable_list
        dictionary[variable] = filter_data(data, variable, class, ["pump_vacuum", "signal_vacuum", "intermediate_vacuum"])
    end
    
    # save and return data
    HarmonicBalance.save(sweep_directory * "/dict" *string(result_num), dictionary)
    return dictionary
    
end 



# =================================================================================================
# function to assemble dictionary for all harmonic variables for all results in a 3-parameter sweep
# =================================================================================================



function create_all_dictionaries(sweep_directory, class)

    # collect a list of all the filenames in the given directory
    file_list = readdir(sweep_directory) 

    # Initialize an empty list
    results_list = []  

    # loop through all files in file_list and append all filenames that begin with "result" to results_list
    for filename in file_list
        if startswith(filename, "results")
            push!(results_list, filename)
        end 
    end 
    
    # collect number of flux points 
    flux_num = size(results_list)[1] 

    # for loop to collect dictionaries for each simulation in flux sweep
    for result_num in 1:(flux_num + 1)
        try
            create_dictionary(sweep_directory, result_num, class)
        catch
            if result_num < (flux_num + 1)
                println("No results found at $sweep_directory/result$result_num.jld2")
            end
            continue
        end
    end 
    
    # return the number of flux points 
    return flux_num
end 