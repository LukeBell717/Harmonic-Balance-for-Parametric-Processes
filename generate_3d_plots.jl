# this script contains functions that help generate only 3D plots

function specify_xy_axes(x_axis, y_axis, f_array, f_cut, E_p_array, E_p_cut, Φ_array, Φ_cut)
    """
    function that specifies array that will be y-axis, and notifies 
    user if there are erroneous cut values or if erroneous value for "y_axis"
    """

    # find array for x-axis
    if x_axis == "frequency"
        if isfinite(f_cut)
            println("Do not specify value for f_cut")
            return
        elseif (isnan(E_p_cut) && isfinite(Φ_cut)) || (isfinite(E_p_cut) && isnan(Φ_cut))
            x = f_array
            xlabel = "Frequency [GHz]"
        else
            println("Please specify one cut value for either E_p_cut and Φ_cut")
            return
        end
    elseif x_axis == "power"
        if isfinite(E_p_cut)
            println("Do not specify value for E_p_cut")
            return
        elseif (isnan(f_cut) && isfinite(Φ_cut)) || (isfinite(f_cut) && isnan(Φ_cut))
            x = E_p_array
            xlabel = "Power [dBm]"
        else
            println("Please specify one cut value for either f_cut and Φ_cut")
            return
        end
    elseif x_axis == "flux"
        if isfinite(Φ_cut)
            println("Do not specify value for Φ_cut")
            return
        elseif (isnan(f_cut) && isfinite(E_p_cut)) || (isfinite(f_cut) && isnan(E_p_cut))
            x = Φ_array
            xlabel = "Flux [h/2e]"
        else
            println("Please specify one cut value for either f_cut and E_p_cut")
            return
        end
    else 
        println("Please enter valid x-axis such as \"frequency,\" \"power,\" or \"flux\"")
    end 

    # collect data and label for y_axis
    if y_axis == x_axis
        println("Please specify valid y_axis that is different than x_axis")
        return
    elseif y_axis == "frequency" 
        if isnan(f_cut) 
            if x_axis == "power" 
                if isnan(Φ_cut)
                    println("Please specify value for Φ_cut")
                    return
                else 
                    y = f_array
                    ylabel = "Frequency [GHz]"
                end
            elseif x_axis == "flux"
                if isnan(E_p_cut)
                    println("Please specify value for E_p_cut")
                    return
                else 
                    y = f_array
                    ylabel = "Frequency [GHz]"
                end
            end
        else
            println("Please do not specify value for f_cut")
            return
        end
    elseif y_axis == "power" 
        if isnan(E_p_cut)
            if x_axis == "frequency"
                if isnan(Φ_cut)
                    println("Please specify value for Φ_cut")
                    return
                else 
                    y = E_p_array
                    ylabel = "Power [dBm]"
                end
            elseif x_axis == "flux"
                if isnan(f_cut)
                    println("Please specify value for f_cut")
                    return
                else
                    y = E_p_array
                    ylabel = "Power [dBm]"
                end
            end
        else 
            println("Please do not specify value for E_p_cut")
            return
        end
    elseif y_axis == "flux"
        if isnan(Φ_cut)
            if x_axis == "frequency"
                if isnan(E_p_cut)
                    println("Please specify value for E_p_cut")
                    return
                else 
                    y = Φ_array
                    ylabel = "Flux [h/2e]"
                end
            elseif x_axis == "power"
                if isnan(f_cut)
                    println("Please specify value for f_cut")
                    return
                else
                    y = E_p_array
                    ylabel = "Power [dBm]"
                end
            end
        else
            println("Please do not specify value for Φ_cut")
            return
        end 
    else 
        println("Please enter valid y-axis such as \"frequency,\" \"power,\" or \"flux.\"")
        return
    end 

    # collect data for title 
    if isfinite(f_cut)
        title = "Frequency = $f_cut [GHz]"
    elseif isfinite(E_p_cut)
        title = "Power = $E_p_cut [dBm]"
    elseif isfinite(Φ_cut)
        title = "Flux = $Φ_cut [h/2e]"
    end

    # return y_axis and corresponding label for y_axis
    return x, xlabel, ylabel, y, title
end

function generate_heatmaps(photon_class_list, x_axis, xlabel, y_axis, ylabel, 
                           pump_mags_array, pump_phases_array, 
                           signal_mags_array, signal_phases_array, 
                           intermediate_mags_array, intermediate_phases_array)

    """
    function that plots heatmaps of the magnitude and phase of pump, signal, and intermediate photons according 
    to the classes specified by "photons"; accepts two swept variables plotted along x and y axes and
    one cut value for the third variable
    """

    # define an empty plot array
    plot_array = []
    
    if (x_axis == "power" && y_axis == "frequency") || (x_axis == "flux" && y_axis == "frequency") || (x_axis == "flux" && y_axis == "power") 
    
        # if the x and y axes are not within the order of the dictionary (branch x frequency x power x flux)
        # then take the transpose of these matrices to send "x to y" and "y to x"
        pump_mags_array = transpose(pump_mags_array)
        pump_phases_array = transpose(pump_phases_array)
        signal_mags_array = transpose(signal_mags_array)
        signal_phases_array = transpose(signal_phases_array)
        intermediate_mags_array = transpose(intermediate_mags_array)
        intermediate_phases_array = transpose(intermediate_phases_array)
    end
    
    # generate heatmaps
    if "pump" in photon_class_list
    
        # execute for loop that overlays all heatmaps (each corresponding to a different branch) for pump magnitude
        heatmap(pump_mags_array[1], color=:coolwarm, title="Pump Magnitude, " * title, 
                xlabel=xlabel, ylabel=ylabel, fontfamily = "Computer Modern")
        for i in 2:length(pump_mags_array)
            heatmap!(pump_mags_array[i], color=:coolwarm, title="Pump Magnitude, " * title, 
                     xlabel=xlabel, ylabel=ylabel, fontfamily = "Computer Modern")
        end
    
        # append p1 to plot_array
        p1 = plot!()
        push!(plot_array, p1) 
    
        # execute for loop that overlays all heatmaps (each corresponding to a different branch) for pump phases
        heatmap(pump_phases_array[1], color=:coolwarm, title="Pump Phase, " * title, 
                xlabel=xlabel, ylabel=ylabel, fontfamily = "Computer Modern")
        for i in 2:length(pump_phases_array)
            heatmap!(pump_phases_array[i], color=:coolwarm, title="Pump Phase, " * title, 
                     xlabel=xlabel, ylabel=ylabel, fontfamily = "Computer Modern")
        end
    
        # append p2 to plot_array
        p2 = plot!()
        push!(plot_array, p2) 
    end
    
    if "signal" in photon_class_list
        
        # execute for loop that overlays all heatmaps (each corresponding to a different branch) for signal magnitude
        heatmap(signal_mags_array[1], color=:coolwarm, title="Signal Magnitude, " * title, 
                xlabel=xlabel, ylabel=ylabel, fontfamily = "Computer Modern")
        for i in 2:length(signal_mags_array)
            heatmap!(signal_mags_array[i], color=:coolwarm, title="Signal Magnitude, " * title, 
                     xlabel=xlabel, ylabel=ylabel, fontfamily = "Computer Modern")
        end
    
        # append p3 to plot_array
        p3 = plot!()
        push!(plot_array, p3)
    
        # create empty arrays in which to store heatmap data for signal photon order parameter (η_signal)
        η_signal_array = []
        
        # loop through all branches in signal_phases_array and extract the order parameter (η_signal)
        for i in 1:length(signal_phases_array)
            
            # collect η for each coordinate pair on heatmap
            η_signal = mod.(3 .* signal_phases_array[i], 2*π) ./3
            push!(η_signal_array, η_signal)
        end
    
        # execute for loop that overlays all heatmaps (each corresponding to a different branch) for signal phase
        heatmap(η_signal_array[1], color=:coolwarm, title="Signal Phase, " * title, 
                xlabel=xlabel, ylabel=ylabel, fontfamily = "Computer Modern")
        for i in 2:length(η_signal_array)
            heatmap!(η_signal_array[i], color=:coolwarm, title="Signal Phase, " * title, 
                     xlabel=xlabel, ylabel=ylabel, fontfamily = "Computer Modern")
        end
    
       # append p4 to plot_array
        p4 = plot!()
        push!(plot_array, p4)
    end
    
    if "intermediate" in photon_class_list
        
        # execute for loop that overlays all heatmaps (each corresponding to a different branch) for intermediate magnitude
        heatmap(intermediate_mags_array[1], color=:coolwarm, title="Intermediate Magnitude, " * title, 
                xlabel=xlabel, ylabel=ylabel, fontfamily = "Computer Modern")
        for i in 2:length(intermediate_mags_array)
            heatmap!(intermediate_mags_array[i], color=:coolwarm, title="Intermediate Magnitude, " * title, 
                     xlabel=xlabel, ylabel=ylabel, fontfamily = "Computer Modern")
        end
    
        # append p3 to plot_array
        p5 = plot!()
        push!(plot_array, p5)
        
        # create empty arrays in which to store heatmap data for signal photon order parameter (η_signal)
        η_intermediate_array = []
        
        # loop through all branches in signal_phases_array and extract the order parameter (η_signal)
        for i in 1:length(intermediate_phases_array)
            
            # collect η for each coordinate pair on heatmap
            η_intermediate = mod.(3 .* intermediate_phases_array[i], 2*π) ./3
            push!(η_intermediate_array, η_intermediate)
        end
    
        # execute for loop that overlays all heatmaps (each corresponding to a different branch) for signal phase
        heatmap(η_intermediate_array[1], color=:coolwarm, title="Intermediate Phase, " * title, 
                xlabel=xlabel, ylabel=ylabel, fontfamily = "Computer Modern")
        for i in 2:length(η_intermediate_array)
            heatmap!(η_intermediate_array[i], color=:coolwarm, title="Intermediate Phase, " * title, 
                     xlabel=xlabel, ylabel=ylabel, fontfamily = "Computer Modern")
        end
    
       # append p4 to plot_array
        p6 = plot!()
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
         plot(plot_array..., layout = (1, 2), size = (1100, 400))
    elseif length(plot_array) == 4
         plot(plot_array..., layout = (2, 2), size = (1000, 700))
    elseif length(plot_array) == 6
         plot(plot_array..., layout = (3, 2), size = (1000, 1100))
    end
end