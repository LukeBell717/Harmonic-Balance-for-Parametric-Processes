# import relevant packages and define variables for our eom
using HarmonicBalance, CSV, DataFrames, DataStructures
@variables t, ϕ(t), Π, E_p, γ, f, f_0, u_2, u_3, u_4, u_5, u_6



# ============================
# function to create directory 
# ============================



function create_directory(main_folder, partition, days, hours, minutes, number_of_threads, number_of_blobs)
    
    # specify the directory for sweep data
    if number_of_blobs == 2
        blob_directory = main_folder*"/bi_sweeps"
        prefix = "bi_sweep"
    elseif number_of_blobs == 3 
        blob_directory = main_folder*"/tri_sweeps"
        prefix = "tri_sweep"
    elseif number_of_blobs == 4
        blob_directory = main_folder*"/tetra_sweeps"
        prefix = "tetra_sweep"
    else 
        println("Please specify a valid number of blobs, i.e., an integer 2, 3, or 4")
    end 
    
    # obtain list of all files in sweep directory and create empty list
    file_list = readdir(blob_directory)

    if isempty(file_list)
        sweep_directory_number = 1
    else 
        index_list = [] 
        # loop through all filenames in file_list
        for filename in file_list
            if startswith(filename, prefix)
                # find string of digits in each filename and convert string to integer; 
                # add this integer to index_list
                string_digits = match(r"\d+", filename).match
                index = parse(Int, string_digits)
                push!(index_list, index)
            end
        end 
        # determine sweep number
        sweep_directory_number = maximum(index_list) + 1
    end 
    
    # convert the integer hours into acceptable string
    if hours < 10
        hours = "0$hours"
    else 
        hours = "$hours"
    end

    # convert the integer minutes into acceptable string
    if minutes < 10 
        minutes = "0$minutes"
    else 
        minutes = "$minutes"
    end

    # specify estimated run time as acceptable string
    if days == 0 
        time = hours * ":" * minutes * ":00"
    elseif days > 0
        time = "$days-" * hours * ":" * minutes * ":00"
    else
        println("Please enter a valid number of days")
    end

    # make new directory in which to store data 
    # copy folder containing execute_sweep_on_grace from script directory to sweep directory
    sweep_directory = blob_directory*"/"*prefix*"$sweep_directory_number"
    mkdir(sweep_directory)
    cp(main_folder*"/scripts/execute_sweep_on_grace.jl", sweep_directory*"/execute_sweep_on_grace.jl")
    
    # generate line to insert into terminal which specifies resouce allocation for Grace
    resource_allocation = "Terminal command to prepare bash script:\n\n" *
                          "cd $sweep_directory; module load dSQ; dsq --job-file executable$sweep_directory_number.txt " *
                          "--time=" * time * " --partition=" * partition * " --nodes=1 --ntasks=1 --cpus-per-task=$number_of_threads \n" 
    println(resource_allocation)
    
    println("Data stored in "*sweep_directory)
    
    return sweep_directory_number, sweep_directory
end



# ==============================================
# function to save sweep parameters in .txt file  
# ==============================================



function save_parameters(main_folder, partition, days, hours, minutes, number_of_threads, number_of_blobs, sweep_directory,
                         sweep_directory_number, outermost_sweep_variable, frequencies, γ_value, u3, u4, u5, u6,
                         f_start, f_end, f_N, E_p_start, E_p_end, E_p_N, Φ_start, Φ_end) 
    
    # collect matrix of expansion coefficients, frequency and flux values for sweep; 
    # note that this must be a matrix with 6 columns 
    nina_matrix = CSV.read(main_folder*"/nina_data/"*"musnail_flux_sweep_u6.csv", DataFrame) 
    
    # convert Φ_start and Φ_end, which are in units of 2*π, into angular units
    Φ_start_converted = 2*π*Φ_start
    Φ_end_converted = 2*π*Φ_end

    # collect the indices that capture the flux between Φ_start and Φ_end
    Φ_data = nina_matrix[:, end] 
    idx0 = argmin(abs.(Φ_data .- Φ_start_converted)) 
    idx1 = argmin(abs.(Φ_data .- Φ_end_converted))
        
    # calculate the total number of flux points in our flux sweep
    Φ_N = length(Φ_data[idx0 : idx1])  
    
    # convert array of frequencies into string interpretable by the terminal
    # for example, the array [f, 3/2*f, 2*f] will be converted to "f, 3/2*f, 2*f"
    frequency_string = join(string.(frequencies), ", ")
    
    # create dictionary to store all input parameters
    parameter_dict = OrderedDict{String, Any}(
        
    # dissipation rate and list of frequencies 
    "sweep_directory_number" => sweep_directory_number,
    "partition" => partition,
    "days" => days, 
    "hours" => hours, 
    "minutes" => minutes,
    "number_of_threads" => number_of_threads, 
    "number_of_blobs" => number_of_blobs,
    "outermost_sweep_variable" => outermost_sweep_variable, 
    "frequencies" => frequency_string * " [f is the signal frequency, measured in GHz]",
    "γ" => "$γ_value" * " [GHz]",  
        
    # expansion coefficients 
    "u3" => u3, 
    "u4" => u4, 
    "u5" => u5, 
    "u6" => u6,    

    # values for frequency sweep
    "f_start" => "$f_start" * " [GHz]",
    "f_end" => "$f_end" * " [GHz]",
    "f_N" => f_N,   
        
    # values for power sweep
    "E_p_start" => "$E_p_start" * " [energy units]",
    "E_p_end" => "$E_p_end" * " [energy units]",
    "E_p_N" => E_p_N,    
        
    # values for flux sweep
    "Φ_start" => "$Φ_start" * " [h/(2*e)]", 
    "Φ_end" => "$Φ_end" * " [h/(2*e)]",
    "Φ_N" => Φ_N
    )
    
    # Open a file in write mode
    parameter_file = open(joinpath(sweep_directory, "parameters$sweep_directory_number.txt"), "w")

    # Write keys and entries to the file
    for (key, value) in parameter_dict
        println(parameter_file, "$key: $value")
    end

    # Close the file
    close(parameter_file)
    
    return Φ_N
    
    # specify parameters keys and keyword arguments
    parameter_keys = ["partition", "days", "hours", "minutes", "number_of_threads", "number_of_blobs",
               "outermost_sweep_variable", "frequencies", "γ_value", "u3", "u4", "u5", "u6",
               "f_start", "f_end", "f_N", "E_p_start", "E_p_end", "E_p_N", "Φ_start", "Φ_end", "Φ_N"]
    parameter_values = [partition, days, hours, minutes, number_of_threads, number_of_blobs,
          outermost_sweep_variable, frequencies, γ_value, u3, u4, u5, u6, 
          f_start, f_end, f_N, E_p_start, E_p_end, E_p_N, Φ_start, Φ_end, Φ_N]
    
    # create and save dictionary 
    parameter_dict = Dict(zip(parameter_keys, parameter_values))
    HarmonicBalance.save(sweep_directory*"/parameter_dict$sweep_directory_number", parameter_dict)
end



# ===================================
# function to find harmonic equations
# ===================================



function collect_harmonic_equations(sweep_directory, sweep_directory_number, number_of_blobs, frequencies, u3, u4, u5, u6)
    
    # define 2*Π
    inv_twopi_sq = 1/((2*π)^2) 
    
    # specify pump term eom depending on bifurcation/trifurcation/tetrafurcation
    if number_of_blobs == 2
        pump_term = (1/u_2)*E_p*cos(2*Π*(2*f/3*t))
    elseif number_of_blobs == 3
        pump_term = (1/u_2)*E_p*cos(2*Π*(3*f/2*t))
    elseif number_of_blobs == 4 
        pump_term = (1/u_2)*E_p*cos(2*Π*(4*f*t))
    end
        
    # determine the desired eom based on the arguments u3, u4, and u5
    if u3 == u4 == u5 == u6 == 1
        ϕ_symbolic = (d(ϕ, t, 2))/((2*Π*f_0)^2) + ((2*Π*γ)*d(ϕ,t))/((2*Π*f_0)^2) + ϕ + pump_term +
                     (u_3/(2*u_2))*ϕ^2 + (u_4/(6*u_2))*ϕ^3 + (u_5/(24*u_2))*ϕ^4 + (u_6/(120*u_2))*ϕ^5 
    elseif u3 == u4 == u5 == 1 && u6 == 0
        ϕ_symbolic = (d(ϕ, t, 2))/((2*Π*f_0)^2) + ((2*Π*γ)*d(ϕ,t))/((2*Π*f_0)^2) + ϕ + pump_term +
                     (u_3/(2*u_2))*ϕ^2 + (u_4/(6*u_2))*ϕ^3 + (u_5/(24*u_2))*ϕ^4 
    elseif u4 == u5 == 1 && u3 == u6 == 0
        ϕ_symbolic = (d(ϕ, t, 2))/((2*Π*f_0)^2) + ((2*Π*γ)*d(ϕ,t))/((2*Π*f_0)^2) + ϕ + pump_term +
                     (u_4/(6*u_2))*ϕ^3 + (u_5/(24*u_2))*ϕ^4
    elseif u3 == u4 == 1 && u5 == u6 == 0 
        ϕ_symbolic = (d(ϕ, t, 2))/((2*Π*f_0)^2) + ((2*Π*γ)*d(ϕ,t))/((2*Π*f_0)^2) + ϕ + pump_term + 
                     (u_3/(2*u_2))*ϕ^2 + (u_4/(6*u_2))*ϕ^3 
    elseif u3 == 1 && u4 == u5 == u6 == 0
        ϕ_symbolic = (d(ϕ, t, 2))/((2*Π*f_0)^2) + ((2*Π*γ)*d(ϕ,t))/((2*Π*f_0)^2) + ϕ + pump_term + 
                     (u_3/(2*u_2))*ϕ^2  
    else
        ϕ_symbolic = "Please enter valid eom"
    end
    
    # display symbolic eom
    display(ϕ_symbolic)
    
    # define as differential equation
    ϕ_eom = DifferentialEquation(ϕ_symbolic, ϕ)
    
    angular_frequencies = 2*Π*frequencies    
    # specify the ansatz ϕ = u(T) cos(2*Π*f*t) + v(T) sin(2*Π*f*t) +...
    add_harmonic!(ϕ_eom, ϕ, angular_frequencies)

    # convert the differential equation to the algebraic harmonic equations
    ϕ_harmonic_equations = get_harmonic_equations(ϕ_eom)
    
    HarmonicBalance.save(sweep_directory*"/harmonic_equations$sweep_directory_number", ϕ_harmonic_equations)
    
    return ϕ_harmonic_equations
end



# ================================================
# function to generate bash script for the cluster
# ================================================



function generate_bash_script(main_folder, number_of_threads, sweep_directory, sweep_directory_number, γ_value,
                              f_start, f_end, f_N, E_p_start, E_p_end, E_p_N, Φ_start, Φ_end, Φ_N) 
    
    # list of parameters that will be passed as command-line arguments     
    parameters = ("\"$main_folder\" " * "\"$sweep_directory\" " * "$sweep_directory_number " * "$γ_value " 
                  * "$f_start " *"$f_end " *"$f_N " * "$E_p_start " * "$E_p_end " * "$E_p_N "
                  * "$Φ_start " *"$Φ_end ") 
    
    command_line_args = ("module load StdEnv; cd $main_folder/scripts; export NUM_JULIA_THREADS=$number_of_threads; " 
                       * "module load Julia; julia execute_sweep_on_grace.jl ")
    
    # open bash script in write mode 
    bash_file = open(joinpath(sweep_directory, "executable$sweep_directory_number.txt"), "w")
    
    # print all command line statements to run sweep over power and frequency for each Φ_N flux points
    for sweep_index in 1:Φ_N
        println(bash_file, command_line_args * "$sweep_index " * parameters)
    end
    close(bash_file)
end



# ===============================================================
# function to prepare all directories and data to run simulations
# ===============================================================



function prepare_sweep_linear(main_folder, partition, days, hours, minutes, number_of_threads, number_of_blobs, outermost_sweep_variable,                                 frequencies, γ_value, u3, u4, u5, u6, f_start, f_end, f_N, E_p_start, E_p_end, E_p_N, Φ_start, Φ_end)          
    
    # create directory in which to store data
    sweep_directory_number, sweep_directory = create_directory(main_folder, partition, days, hours,
                                                               minutes, number_of_threads, number_of_blobs)
    
    # save data into text file in corresponding directory
    Φ_N = save_parameters(main_folder, partition, days, hours, minutes, number_of_threads, number_of_blobs, sweep_directory,
                          sweep_directory_number, outermost_sweep_variable, frequencies, γ_value, u3, u4, u5, u6,
                          f_start, f_end, f_N, E_p_start, E_p_end, E_p_N, Φ_start, Φ_end)
    
    # collect harmonic equations for sweep
    ϕ_harmonic_equations = collect_harmonic_equations(sweep_directory, sweep_directory_number, 
                                                      number_of_blobs, frequencies, u3, u4, u5, u6)
    
    # generate bash script to run on Grace 
    generate_bash_script(main_folder, number_of_threads, sweep_directory, sweep_directory_number, γ_value,
                         f_start, f_end, f_N, E_p_start, E_p_end, E_p_N, Φ_start, Φ_end, Φ_N) 
    
    return ϕ_harmonic_equations
end