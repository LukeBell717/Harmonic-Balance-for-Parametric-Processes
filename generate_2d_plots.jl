# this script contains functions that generate both 2D and/or 3D plots

function specify_x_axis(x_axis, f_array, f_cut, E_p_array, E_p_cut, Φ_array, Φ_cut)
    """
    function that specifies array that will be x-axis, and notifies 
    user if there are erroneous cut values 
    """
    
    if x_axis == "frequency"
        if any(isnan, [E_p_cut, Φ_cut])
            println("Please specify values for E_p_cut and Φ_cut")
            return
        elseif isfinite(f_cut)
            println("Do not specify value for f_cut")
            return
        else
            x = f_array
            xlabel = "Frequency [GHz]"
        end
    elseif x_axis == "power"
        if any(isnan, [f_cut, Φ_cut])
            println("Please specify values for f_cut and Φ_cut")
            return
        elseif isfinite(E_p_cut)
            println("Do not specify value for E_p_cut")
            return
        else
            x = E_p_array
            xlabel = "Power [dBm]"
        end
    elseif x_axis == "flux"
        if any(isnan, [f_cut, E_p_cut])
            println("Please specify values for f_cut and E_p_cut")
            return
        elseif isfinite(Φ_cut)
            println("Do not specify value for Φ_cut")
            return
        else
            x = Φ_array
            xlabel = "Flux [h/2e]"
        end
    else 
        println("Please enter valid x-axis such as \"frequency,\" \"power,\" or \"flux.\"")
    end 

    return x, xlabel
end

function find_cut_indices(f_array, f_cut, E_p_array, E_p_cut, Φ_array, Φ_cut)
    
    """ 
    the default for f_cut, E_p_cut and Φ_cut is NaN; 
    if we specify a non-NaN value for one of these cuts then
    we find the index of this cut in the corresponding array
    """

    # if necessary, find index at which to cut f_array
    if isnan(f_cut)
        f_cut_index = Symbol(":")
    elseif isfinite(f_cut)
        f_cut_index = argmin(abs.(f_array .- f_cut))
    else 
        println("Please specify valid f_cut")
    end 

    # if necessary, find index at which to cut E_p_array
    if isnan(E_p_cut)
        E_p_cut_index = Symbol(":")
    elseif isfinite(E_p_cut)
        E_p_cut_index = argmin(abs.(E_p_array .- E_p_cut))
    else 
        println("Please specify valid E_p_cut")
    end

    # if necessary, find index at which to cut Φ_array; 
    # note that we multiply Φ_cut by 2*π since we must convert from
    # angular to linear frequency
    if isnan(Φ_cut)
        Φ_cut_index = Symbol(":")
    elseif isfinite(Φ_cut)
        Φ_cut_index = argmin(abs.(Φ_array .- 2*π*Φ_cut))
    else 
        println("Please specify valid Φ_cut")
    end
    
    return f_cut_index, E_p_cut_index, Φ_cut_index
end


function collect_plot_data(photon_class_list, full_dict, branch_N, f_cut_index, E_p_cut_index, Φ_cut_index)
    
    """ 
    collects data for the magnitude and phase of pump, signal, and intermediate photons according 
    to the classes specified by "photons"
    """

    # define empty arrays in which to store solutions
    # for magnitudes and phases of pump, signal, and intermediate photons
    # the length of array corresponds to the number of branches 
    pump_mags_array= []
    pump_phases_array = []
    signal_mags_array = []
    signal_phases_array = []
    intermediate_mags_array = []
    intermediate_phases_array = []
    
    # define empty arrays to store harmonic coefficients
    # u1_array = []
    # v1_array = []
    # u2_array = []
    # v2_array = []
    # u3_array = []
    # v3_array = []

    # loops through all branches and computes magnitude and phase of pump, signal, and intermediate 
    # photons for indices of branch, frequency, power, and flux points
    for branch_index in 1:branch_N

        # compute magnitudes and phases of pump photons if we specify "pump" in our list of photons
        if "pump" in photon_class_list

            # collect harmonic coefficients of the pump photons
            u1 = full_dict["u1"][branch_index, eval(f_cut_index), eval(E_p_cut_index), eval(Φ_cut_index)]
            v1 = full_dict["v1"][branch_index, eval(f_cut_index), eval(E_p_cut_index), eval(Φ_cut_index)]

            # find magnitude and phase of pump photons
            pump_mag = sqrt.(u1.^2 + v1.^2)
            pump_phase = atan.(v1, u1) 

            # append pump_mag (either array or matrix) to pump_mags_array and 
            # pump_phase (either array or matrix) to pump_phases_array if 
            # both pump_mag and pump_phase contain at least one non-NaN value
            if any(isfinite, pump_mag) && any(isfinite, pump_phase)
                push!(pump_mags_array, pump_mag)
                push!(pump_phases_array, pump_phase)
                # push!(u1_array, u1)
                # push!(v1_array, v1)
            end 
        end

        # compute magnitudes and phases of signal photons if we specify "signal" in our list of photons
        if "signal" in photon_class_list

            # collect harmonic coefficients of the signal photons
            u2 = full_dict["u2"][branch_index, eval(f_cut_index), eval(E_p_cut_index), eval(Φ_cut_index)]
            v2 = full_dict["v2"][branch_index, eval(f_cut_index), eval(E_p_cut_index), eval(Φ_cut_index)]

            # compute magnitude and phase of signal photons
            signal_mag = sqrt.(u2.^2 + v2.^2)
            signal_phase = atan.(v2, u2) 

            # append signal_mag (either array or matrix) to signal_mags_array and 
            # signal_phases (either array or matrix) to signal_phases_array if 
            # both signal_mag and signal_phase contain at least one non-NaN value
            if any(isfinite, signal_mag) && any(isfinite, signal_phase)
                push!(signal_mags_array, signal_mag)
                push!(signal_phases_array, signal_phase)
                # push!(u2_array, u2)
                # push!(v2_array, v2)
            end 
        end

        # compute magnitudes and phases of intermediate photons if we specify "intermediate" in our list of photons
        if "intermediate" in photon_class_list

            # collect harmonic coefficients of the intermediate photons
            u3 = full_dict["u3"][branch_index, eval(f_cut_index), eval(E_p_cut_index), eval(Φ_cut_index)]
            v3 = full_dict["v3"][branch_index, eval(f_cut_index), eval(E_p_cut_index), eval(Φ_cut_index)]

            # compute magnitude and phase of intermediate photons
            intermediate_mag = sqrt.(u3.^2 + v3.^2)
            intermediate_phase = atan.(v3, u3) 

            # append intermediate_mag (either array or matrix) to intermediate_mags_array 
            # and intermediate_phases (either array or matrix) to intermediate_phases_array if 
            # both intermediate_mag and intermediate_phase contain at least one non-NaN value
            if any(isfinite, intermediate_mag) && any(isfinite, intermediate_phase)
                push!(intermediate_mags_array, intermediate_mag)
                push!(intermediate_phases_array, intermediate_phase)
                # push!(u3_array, u3)
                # push!(v3_array, v3)
            end
        end
    end
    
    return pump_mags_array, pump_phases_array, signal_mags_array, 
           signal_phases_array, intermediate_mags_array, intermediate_phases_array 
           # u1_array, v1_array, u2_array, v2_array, u3_array, v3_array
end



function generate_plots(photon_class_list, x, xlabel, 
                   pump_mags_array, pump_phases_array, 
                   signal_mags_array, signal_phases_array, 
                   intermediate_mags_array, intermediate_phases_array) 

    """
    function that plots the magnitude and phase of pump, signal, and intermediate photons according 
    to the classes specified by "photons"
    """

    # define an empty plot array
    plot_array = []

    # plot data and specify graph characteristics
    if "pump" in photon_class_list
        p1 = plot(x, pump_mags_array, xlabel = xlabel, ylabel = "ϕ", title = "Pump Magnitude", fontfamily="Computer Modern")
        # plot!(legend=:topright)
        p2 = plot(x, pump_phases_array, xlabel = xlabel, ylabel = "rad", title = "Pump Phase", fontfamily="Computer Modern")
        # plot!(legend=:topright)
        push!(plot_array, p1)
        push!(plot_array, p2)
    end
    
    if "signal" in photon_class_list
        p3 = plot(x, signal_mags_array, xlabel = xlabel , ylabel = "ϕ", title = "Signal Magnitude", fontfamily="Computer Modern")
        # plot!(legend=:topleft)
        p4 = plot(x, signal_phases_array, xlabel = xlabel , ylabel = "rad", title = "Signal Phase", fontfamily="Computer Modern")
        # plot!(legend=:topleft)
        push!(plot_array, p3)
        push!(plot_array, p4)
    end
    
    if "intermediate" in photon_class_list
        p5 = plot(x, intermediate_mags_array, xlabel = xlabel , ylabel = "ϕ", title = "Intermediate Magnitude", fontfamily="Computer Modern")
        #plot!(legend=:topleft)
        p6 = plot(x, intermediate_phases_array, xlabel = xlabel , ylabel = "rad", title = "Intermediate Phase", fontfamily="Computer Modern")
        #plot!(legend=:topleft)
        push!(plot_array, p5)
        push!(plot_array, p6)
    end

    # create empty list into which we will store the number of branches in each plot
    n_values = Int[]
    
    for plot in plot_array 
        # convert name of plot into a string
        plot_string = string(plot)
    
        # search for a pattern in a string where it expects 
        # to find "n=" followed by one or more digits
        regex_pattern = r"n=(\d+)"
    
        # use the match function to find the number in the string
        match_result = match(regex_pattern, plot_string)
    
        if match_result !== nothing
            # Extract the matched number as a string
            branches_number = push!(n_values, parse(Int, match_result.captures[1]))
        end
    end

    # notify if there are no solutions
    if all(x -> x == 1, n_values)
        println("There are no solutions for this sweep over " * x_axis * " given your cut values.")
    end

    # plot data according to the classes of photons specified as argument
    if length(plot_array) == 2
         plot(plot_array..., layout = (1, 2), size = (950, 400))
    elseif length(plot_array) == 4
         plot(plot_array..., layout = (2, 2), size = (950, 700))
    elseif length(plot_array) == 6
         plot(plot_array..., layout = (3, 2), size = (900, 1000))
    end
end