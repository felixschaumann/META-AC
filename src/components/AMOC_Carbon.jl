@defcomp AMOC_Carbon begin

    # Variables
    AMOC_decrease = Variable(index=[time], unit="Sv") # AMOC decrease (Sv) - a positive value indicates a decrease in AMOC strength
    AMOC_strength = Variable(index=[time], unit="Sv") # AMOC strength (Sv)
    CO2_AMOC = Variable(index=[time], unit="GtCO2") # CO2 release due to AMOC decrease (GtCO2)
    cum_CO2_AMOC = Variable(index=[time], unit="GtCO2")
    
    # Parameters
    T_AT = Parameter(index=[time], unit="degC")
    std_AMOC_strength = Parameter(index=[time], default=zeros(451), unit="Sv") # standard deviation of AMOC strength in a given projection
    
    beta_AMOC = Parameter(index=[time]) # coefficient representing the sensitivity of AMOC to temperature
    eta_AMOC = Parameter() # coefficient representing the sensitivity of the carbon drawdown to AMOC decrease
    
    T_AT_start = Parameter(default=0.854, unit="degC")
    AMOC_start = Parameter(default=19.0, unit="Sv")
    MPI_AMOC_pi = Parameter(unit="Sv") # AMOC strength in MPI-ESM1.2-LR in preindustrial simulation
    AMOC_pi = Parameter(unit="Sv", default=20.) # preindustrial AMOC strength of respective model
    std_AMOC_pi = Parameter(unit="Sv", default=0.) # standard deviation of preindustrial AMOC strength

    function run_timestep(pp, vv, dd, tt)
        
        # AMOC decrease (measured in Sv)
        vv.AMOC_decrease[tt] = pp.beta_AMOC[tt] * (pp.T_AT[tt] - pp.T_AT_start)
        vv.AMOC_strength[tt] = pp.AMOC_start - vv.AMOC_decrease[tt] + pp.std_AMOC_strength[tt]
        
        # calculate CO2 release depending on AMOC strength while not imposing a maximum CO2 release
        vv.CO2_AMOC[tt] = pp.eta_AMOC * (pp.AMOC_pi + pp.std_AMOC_pi - vv.AMOC_strength[tt]) / pp.AMOC_pi * pp.MPI_AMOC_pi
        
        if is_first(tt)    
            vv.cum_CO2_AMOC[tt] = vv.CO2_AMOC[tt]
        else
            vv.cum_CO2_AMOC[tt] = vv.cum_CO2_AMOC[tt-1] + vv.CO2_AMOC[tt]
        end
    end
end

function addAMOC_Carbon(model::Model, calibration::Union{Float64, Int, String}; before=nothing, after=:PatternScaling, scenario::String="ssp245", T_AT_2100=2.289, AMOC_2010=19.0)

    amoc_carbon = add_comp!(model, AMOC_Carbon, before=before, after=after, first=2010)
    
    # calibrate stylised AMOC decrease from 2010 to 2100
    if typeof(calibration) == Float64 || typeof(calibration) == Int
        T_AT_2010 = 0.854
        amoc_carbon[:AMOC_start] = AMOC_2010
        println("Calibrating AMOC to stylised 2100 decrease of $calibration Sv, starting from an AMOC strength of $(AMOC_2010) Sv in 2010.")
        amoc_carbon[:beta_AMOC] = (calibration / (T_AT_2100 - T_AT_2010)) * ones(dim_count(model, :time))
        amoc_carbon[:MPI_AMOC_pi] = 18.9558
        amoc_carbon[:AMOC_pi] = AMOC_2010 # assume that AMOC strength in 2010 is the same as in preindustrial times
        amoc_carbon[:std_AMOC_pi] = 0.0
        amoc_carbon[:std_AMOC_strength] = zeros(dim_count(model, :time))
    # calibrate to projection from 2015 to 2100
    elseif typeof(calibration) == String
        println("Calibrating AMOC to follow $calibration projection for $(scenario) from 2015 to 2100.")
        beta_AMOC, AMOC_2015, T_AT_2015, AMOC_pi, MPI_AMOC_pi = calibrate_AMOC(model, calibration; scenario=scenario)
        amoc_carbon[:beta_AMOC] = beta_AMOC
        amoc_carbon[:AMOC_start] = AMOC_2015
        amoc_carbon[:T_AT_start] = T_AT_2015
        amoc_carbon[:AMOC_pi] = AMOC_pi
        amoc_carbon[:MPI_AMOC_pi] = MPI_AMOC_pi
        amoc_carbon[:std_AMOC_pi] = 0.0 # will be updated during Monte Carlo runs
        amoc_carbon[:std_AMOC_strength] = zeros(dim_count(model, :time)) # will be updated during Monte Carlo runs
    end
    
    # calibrate carbon cycle effect
    eta_AMOC = 0.023 * 22/6 # is updated by datadir() as output of hosing_regressions.jl when calling with get_amoc_carbon_model()
    amoc_carbon[:eta_AMOC] = eta_AMOC # will be updated during Monte Carlo runs

    amoc_carbon
end