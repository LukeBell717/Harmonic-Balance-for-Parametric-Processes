function analyze_data(sweep_directory, photon_class_list, x_axis, y_axis, f_cut, E_p_cut, Φ_cut)

    """
    function to analyze data for a given sweep by plotting magnitude and phase of 
    pump, signal, and/or intermediate photons as a function of frequency, power, and/or flux

    Args:
        sweep_directory: string specifying path to directory containing sweep
        photon_class_list: list containing any of the strings "pump", "signal" and/or "intermediate
                            corresponding to photons graphs that will be plotted
        f_cut: where we may cut our sweep over frequency; default value is NaN
        E_p_cut: where we may cut our sweep over power; default value is NaN
        Φ_cut: where we may cut our sweep over flux; default value is NaN

    Returns:
        returns 2, 4, or 6 graphs that are either 2D, corresponding to two cut values 
        or heatmaps, corresponding to one cut value
    """
    
    # collect sweep_directory_number
    sweep_directory_number_str = match(r"\d+$", sweep_directory).match
    sweep_directory_number = parse(Int, sweep_directory_number_str)

    # collect full dictionary of data and dictionary containing sweep parameters
    full_dict = HarmonicBalance.load(sweep_directory * "/full_dict$sweep_directory_number.jld2")
    parameter_dict = HarmonicBalance.load(sweep_directory*"/parameter_dict$sweep_directory_number.jld2")

    # collect size of dictionary value for given key
    branch_N, f_N, E_p_N, Φ_N = size(full_dict["u1"])

    # collect sweep parameters 
    f_start = parameter_dict["f_start"]
    f_end = parameter_dict["f_end"]
    f_N = parameter_dict["f_N"]
    E_p_start = parameter_dict["E_p_start"]
    E_p_end = parameter_dict["E_p_end"]
    E_N = parameter_dict["E_p_N"]
    Φ_start = parameter_dict["Φ_start"]
    Φ_end = parameter_dict["Φ_end"]
    Φ_N = parameter_dict["Φ_N"]

    # define arrays for frequency and power
    f_array = LinRange(f_start, f_end, f_N)
    E_p_array = LinRange(E_p_start, E_p_end, E_p_N)

    # collect array for flux
    nina_matrix = CSV.read(main_folder*"/nina_data/"*"musnail_flux_sweep_u6.csv", DataFrame) 
    Φ_data = nina_matrix[:, end] 
    idx0 = argmin(abs.(Φ_data .- 2*π*Φ_start)) 
    idx1 = argmin(abs.(Φ_data .- 2*π*Φ_end))
    Φ_array  = Φ_data[idx0 : idx1]

    if y_axis == "photons"

        # specify and collect data for x_axis of plot
        x, xlabel = specify_x_axis(x_axis, f_array, f_cut, E_p_array, E_p_cut, Φ_array, Φ_cut)
        
        # find indices over which we should analyze data
        f_cut_index, E_p_cut_index, Φ_cut_index = find_cut_indices(f_array, f_cut, E_p_array, E_p_cut, Φ_array, Φ_cut)
        
        # collect data specified by f_cut_index, E_p_cut_index, Φ_cut_index
        pump_mags_array, pump_phases_array, 
        signal_mags_array, signal_phases_array, 
        intermediate_mags_array, intermediate_phases_array = collect_plot_data(photon_class_list, full_dict, branch_N, 
                                                                               f_cut_index, E_p_cut_index, Φ_cut_index)
        # plot simulation data
        generate_plots(photon_class_list, x, xlabel, 
                       pump_mags_array, pump_phases_array, 
                       signal_mags_array, signal_phases_array, 
                       intermediate_mags_array, intermediate_phases_array)
        
    elseif y_axis == "frequency" || y_axis == "power" || y_axis == "flux"
    
        # specify and collect data for x_axis and y_axis of plot
        x, xlabel, ylabel, y, title = specify_xy_axes(x_axis, y_axis, f_array, f_cut, E_p_array, E_p_cut, Φ_array, Φ_cut)
    
        # find indices over which we should analyze data
        f_cut_index, E_p_cut_index, Φ_cut_index = find_cut_indices(f_array, f_cut, E_p_array, E_p_cut, Φ_array, Φ_cut)
    
        # collect data specified by f_cut_index, E_p_cut_index, Φ_cut_index
        pump_mags_array, pump_phases_array, 
        signal_mags_array, signal_phases_array, 
        intermediate_mags_array, intermediate_phases_array = collect_plot_data(photon_class_list, full_dict, branch_N, 
                                                                               f_cut_index, E_p_cut_index, Φ_cut_index)
    
        # generate heatmaps of simulation data
        generate_heatmaps(photon_class_list, x_axis, xlabel, y_axis, ylabel, 
                      pump_mags_array, pump_phases_array, 
                      signal_mags_array, signal_phases_array, 
                      intermediate_mags_array, intermediate_phases_array)
        
    else 
        println("Please enter valid x-axis such as \"frequency,\" \"power,\" or \"flux\"")
    end
end