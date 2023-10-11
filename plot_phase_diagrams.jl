function plot_phase_diagrams(main_folder, sweep_directory, class, flux_num, Φ_start, Φ_end)
    
    # collect flux data from nina matrix 
    nina_matrix = CSV.read(main_folder*"/nina_data/musnail_flux_sweep_u6.csv", DataFrame) 
    Φ_data = nina_matrix[:, end]; 
    
    # convert Φ_start and Φ_end, which is in units of 2*π, into angular units
    Φ_start_converted = 2*π*Φ_start
    Φ_end_converted = 2*π*Φ_end

    # collect the indices that capture the flux between Φ_start and Φ_end
    Φ_data = nina_matrix[:, end] 
    idx0 = argmin(abs.(Φ_data .- Φ_start_converted)) 
    idx1 = argmin(abs.(Φ_data .- Φ_end_converted))
    Φ_range = Φ_data[idx0: idx1]
    
    # loop through all flux values 
    for sweep_index in 1:flux_num
        # collect flux value
        flux_value = Φ_range[sweep_index]/(2*π)
        try 
            # display results
            results = HarmonicBalance.load(sweep_directory * "/results$sweep_index.jld2")
            println("Plot $sweep_index: Φ => $flux_value")
            display(plot_phase_diagram(results, class = ["stable", "physical"]))
            savefig("$sweep_directory/phase_diagram$sweep_index.png")
        catch
            println("There are no results for Plot $sweep_index where Φ = $flux_value \n")
            continue 
        end 
    end
end 