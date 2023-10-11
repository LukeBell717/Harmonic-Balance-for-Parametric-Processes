using HarmonicBalance, CSV, DataFrames, ArgParse
@variables ω, E_p, Π, κ, u_2, u_3, u_4, u_5, u_6, ω_0



# ===============================================================
# function to execute single flux sweep given necessary arguments  
# ===============================================================



# this function will be called form the terminal; here sweep_index is the nth sweep out of all Φ_N flux points
function execute_sweep_on_notebook(sweep_index, main_folder, sweep_directory, sweep_directory_number, γ_value, 
                                   f_start, f_end, f_N, E_p_start, E_p_end, E_p_N, Φ_start, Φ_end)
    
    # collect matrix of expansion coefficients, frequency and flux values for sweep; 
    # note that this must be a matrix with 6 columns 
    nina_matrix = CSV.read(main_folder*"/nina_data/"*"musnail_flux_sweep_u6.csv", DataFrame) 
    
    # convert γ, f, and Φ, which are in units of 2*π, into angular units
    Φ_start_converted = 2*π*Φ_start
    Φ_end_converted = 2*π*Φ_end

    # collect the indices that capture the flux between Φ_start and Φ_end
    Φ_data = nina_matrix[:, end] 
    idx0 = argmin(abs.(Φ_data .- Φ_start_converted)) 
    idx1 = argmin(abs.(Φ_data .- Φ_end_converted))
    
    # collect the values of u_n and ω0 within the range specified by Φ_start and Φ_end
    submatrix = nina_matrix[idx0 : idx1, 1: end-1]
    
    # collect values of u_n and ω0
    u2_val = submatrix[sweep_index, 1]
    u3_val = submatrix[sweep_index, 2]
    u4_val = submatrix[sweep_index, 3]
    u5_val = submatrix[sweep_index, 4]
    u6_val = submatrix[sweep_index, 5]
    f0_val = submatrix[sweep_index, 6] # note submatrix[sweep_index, 6] is in GHz, so we must convert to angular frequency 

    # define varied and fixed parameters for flux sweep and load harmonic equations from sweep directory
    varied = (f => LinRange(f_start, f_end, f_N), E_p => LinRange(E_p_start, E_p_end, E_p_N))
    fixed = (Π => π, γ => γ_value, u_2 => u2_val, u_3 => u3_val, u_4 => u4_val, u_5 => u5_val, u_6 => u6_val, f_0 => f0_val)
    ϕ_harmonic_equations = HarmonicBalance.load(sweep_directory*"/harmonic_equations$sweep_directory_number"*".jld2")

    try
        # solve ϕ_eom using the varied and fixed parameters; threading = true enables parallelization
        results = get_steady_states(ϕ_harmonic_equations, varied, fixed, threading = true)
        # save data with a filename containing sweep_index
        HarmonicBalance.save(sweep_directory*"/results$sweep_index", results)
    catch 
        # if one of the sweeps throws an error (which has happened at Φ = -π) 
        # then disregard and continue to subsequent sweep
        println("We could not process results$sweep_index")
    end
end