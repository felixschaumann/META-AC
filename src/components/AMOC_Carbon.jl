@defcomp AMOC_Carbon begin

    # Variables
    AMOC_decrease = Variable(index=[time], unit="Sv") # AMOC decrease (Sv) - a positive value indicates a decrease in AMOC strength
    AMOC_strength = Variable(index=[time], unit="Sv") # AMOC strength (Sv)
    CO2_AMOC = Variable(index=[time], unit="GtCO2") # CO2 release due to AMOC decrease (GtCO2)
    cum_CO2_AMOC = Variable(index=[time], unit="GtCO2")
    
    # Parameters
    T_AT = Parameter(index=[time], unit="degC")
    std_AMOC_strength = Parameter(index=[time], default=zeros(451), unit="Sv") # standard deviation of AMOC strength in a given projection
    
    I_AMOC = Variable{Bool}(index=[time])
    deltaT_country_AMOC = Variable(index=[time, country], unit="degC")
    T_country_AMOC = Variable(index=[time, country], unit="degC")
    max_deltaT_country_AMOC = Parameter(index=[country], unit="degC")
    Delta_AMOC = Parameter(unit="year", default=5)
    T_country_base = Parameter(index=[time, country], unit="degC")
    pattern_onset = Parameter(default=2015)

    beta_AMOC = Parameter(index=[time]) # coefficient representing the sensitivity of AMOC to temperature
    eta_AMOC = Parameter() # coefficient representing the sensitivity of the carbon drawdown to AMOC decrease
    
    T_AT_start = Parameter(default=0.854, unit="degC")
    AMOC_start = Parameter(default=19.0, unit="Sv")
    MPI_AMOC_pi = Parameter(unit="Sv") # AMOC strength in MPI-ESM1.2-LR in preindustrial simulation
    AMOC_pi = Parameter(unit="Sv", default=20.) # preindustrial AMOC strength of respective model
    std_AMOC_pi = Parameter(unit="Sv", default=0.) # standard deviation of preindustrial AMOC strength

    function run_timestep(pp, vv, dd, tt)
        
        vv.I_AMOC[tt] = false
        if tt >= TimestepValue(pp.pattern_onset)
            vv.I_AMOC[tt] = true
        end
        
        # AMOC decrease (measured in Sv)
        vv.AMOC_decrease[tt] = pp.beta_AMOC[tt] * (pp.T_AT[tt] - pp.T_AT_start)
        vv.AMOC_strength[tt] = pp.AMOC_start - vv.AMOC_decrease[tt] + pp.std_AMOC_strength[tt]
        
        # calculate CO2 release depending on AMOC strength while not imposing a maximum CO2 release
        vv.CO2_AMOC[tt] = pp.eta_AMOC * (pp.AMOC_pi + pp.std_AMOC_pi - vv.AMOC_strength[tt]) / pp.AMOC_pi * pp.MPI_AMOC_pi
        
        if is_first(tt)    
            vv.cum_CO2_AMOC[tt] = vv.CO2_AMOC[tt]
            for cc in dd.country
                vv.deltaT_country_AMOC[tt, cc] = vv.I_AMOC[tt] ? pp.max_deltaT_country_AMOC[cc] / pp.Delta_AMOC : 0
            end
        else
            vv.cum_CO2_AMOC[tt] = vv.cum_CO2_AMOC[tt-1] + vv.CO2_AMOC[tt]
    
            for cc in dd.country
                if vv.I_AMOC[tt] && abs(vv.deltaT_country_AMOC[tt-1, cc]) < abs(pp.max_deltaT_country_AMOC[cc])
                    vv.deltaT_country_AMOC[tt, cc] = vv.deltaT_country_AMOC[tt-1, cc] + pp.max_deltaT_country_AMOC[cc] / pp.Delta_AMOC
                    if abs(vv.deltaT_country_AMOC[tt, cc]) > abs(pp.max_deltaT_country_AMOC[cc])
                        vv.deltaT_country_AMOC[tt, cc] = pp.max_deltaT_country_AMOC[cc]
                    end
                else
                    vv.deltaT_country_AMOC[tt, cc] = vv.deltaT_country_AMOC[tt-1, cc]
                end
            end
        end

        for cc in dd.country
            vv.T_country_AMOC[tt, cc] = pp.T_country_base[tt, cc] + vv.deltaT_country_AMOC[tt, cc]
        end
    end
end

function addAMOC_Carbon(model::Model, calibration::Union{Float64, Int, String}, temp_pattern::Union{String, Bool}; before=nothing, after=:PatternScaling, scenario::String="ssp245", T_AT_2100=2.289, AMOC_2010=19.0)

    amoc_carbon = add_comp!(model, AMOC_Carbon, before=before, after=after, first=2010)

    if temp_pattern == false
        amoc_carbon[:pattern_onset] = 2500 # set onset year of temperature pattern to beyond model time horizon
        temp_pattern = "IPSL" # add dummy value to avoid error in loading dataset
    end
    
    params = CSV.read("../data/AMOCparams.csv", DataFrame)
    amoc_carbon[:max_deltaT_country_AMOC] = [(iso ∈ params."Country code" ? params[params."Country code" .== iso, temp_pattern][1] : 0.0) for iso in dim_keys(model, :country)]

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