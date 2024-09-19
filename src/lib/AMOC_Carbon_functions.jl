using CSV, DataFrames, NCDatasets, Mimi

include("../components/AMOC_Carbon.jl")
include("../components/CO2Converter_AMOC.jl")

cmip6_model_names = Dict(
    "ACCESS-ESM1-5" => "accessesm15",
    "CanESM5" => "canesm5",
    "CESM2" => "cesm2",
    "MIROC-ES2L" => "miroces2l",
    "MPI-ESM1.2-LR" => "mpiesmlr",
    "NorESM2-LM" => "noresm2lm",
    "GISS-E2-1-G" => "gisse21g",
)

#%%
# FUNCTIONS FOR BUILDING AND RUNNING THE AMOC-CARBON MODEL
##########################################################################################

function get_primary_model_name(model_name::String, model_dict::Dict{String, String})
    for (key, value) in model_dict
        if model_name == key || model_name == value
            return key
        end
    end
    error("Model name not found in CMIP6 dictionary.")
end

function get_projection_ds(model_name::String; scenario::String="ssp245")
    model_name = get_primary_model_name(model_name, cmip6_model_names)    
    model_dir = "../data/CMIP6_amoc/$(model_name)/"
    model_dir_files = readdir(model_dir)
    
    search_string = "msftmz_$(cmip6_model_names[model_name])_$(scenario)"
    scenario_files = filter(x -> occursin(search_string, x), model_dir_files)
    
    if length(model_dir_files) > 1
        try # to allow for single NorESM ensemble members for ssp126 and ssp858
            full_ens_path = filter(x -> !occursin(search_string*"_r1-3i", x), scenario_files)[1]
        catch
            full_ens_path = filter(x -> !occursin(search_string*"_r1i", x), scenario_files)[1]
        end
    else
        full_ens_path = scenario_files[1]
    end
    
    full_amoc_ds = NCDatasets.Dataset(joinpath(model_dir, full_ens_path))
    return full_amoc_ds
end

function calibrate_AMOC(model::Model, esm_calibration::String; scenario::String="ssp245")

    esm_calibration = get_primary_model_name(esm_calibration, cmip6_model_names)
    amoc_ds = get_projection_ds(esm_calibration; scenario=scenario)
    
    global amoc_esm_mean
    try # to allow for single NorESM ensemble members for ssp126 and ssp858
        amoc_esm_mean = amoc_ds["mean"] |> Array
    catch
        amoc_esm_mean = amoc_ds["r1i1p1f1"] |> Array
    end
    
    AMOC_pi_vals = CSV.read("../data/CMIP6_amoc/AMOC_pi_values.csv", DataFrame)
    AMOC_pi = AMOC_pi_vals[:, esm_calibration][1]
    MPI_AMOC_pi = AMOC_pi_vals[:, "MPI-ESM1.2-LR"][1]

    tat = model[:TemperatureConverter, :T_AT][266:351]

    beta_AMOC_2016_2100 = - ((amoc_esm_mean .- amoc_esm_mean[1]) ./ (tat .- tat[1]))[2:end]
    beta_AMOC = vcat(zeros(266) .* beta_AMOC_2016_2100[end], beta_AMOC_2016_2100, ones(100) .* mean(beta_AMOC_2016_2100[end-10:end]))

    return beta_AMOC, amoc_esm_mean[1], tat[1], AMOC_pi, MPI_AMOC_pi
end

function add_AMOC(model, amoc)
    
    amocmodel = addAMOC(model, amoc, after=:PatternScaling);
    
    connect_param!(model, :AMOC=>:T_AT, :TemperatureConverter=>:T_AT);
    connect_param!(model, :AMOC=>:T_country_base, :PatternScaling=>:T_country);
    
    ## Use AMOC temperatures rather than PatternScaling temperatures
    connect_param!(model, :Consumption=>:T_country, :AMOC=>:T_country_AMOC);
    amocmodel[:uniforms] = rand(Uniform(0, 1), dim_count(model, :time));
    
    return model
end

function get_amoc_carbon_model(amoc_calibration; scenario::String="ssp245", Dam::String="pointestimate")
    
    model = base_model(rcp=scenario=="ssp126" ? "RCP3-PD/2.6" : "RCP$(scenario[5]).$(scenario[6])", ssp="SSP$(scenario[4])", tdamage=Dam=="BHM" ? "pointestimate" : Dam, slrdamage="mode"); run(model);
    replace!(model, :CO2Converter => CO2Converter_AMOC);
    Mimi.set_first_last!(model, :CO2Converter, first = 2010, last = 2200);
    set_param!(model, :CO2Converter, :co2_amoc, zeros(dim_count(model, :time)))
    
    addAMOC_Carbon(model, amoc_calibration; scenario=scenario);
    
    connect_param!(model, :AMOC_Carbon=>:T_AT, :TemperatureConverter=>:T_AT);
    connect_param!(model, :CO2Converter=>:co2_amoc, :AMOC_Carbon=>:CO2_AMOC);
    
    update_param!(model, :Consumption, :damagepersist, 0.25);
    run(model)
    
    return model
end

function get_base_model_w_params(d::Dict)
    T_AT_2010, T_AT_2100 = 0.854, 2.289;
    @unpack Dam, AMOC, prtp, persist, Scenario = d
    if AMOC in ["T", "T&C"]
        @unpack Calib = d
    end
    if AMOC=="none"
        m = base_model(rcp=Scenario=="ssp126" ? "RCP3-PD/2.6" : "RCP$(Scenario[5]).$(Scenario[6])", ssp="SSP$(Scenario[4])", tdamage=Dam=="BHM" ? "pointestimate" : Dam, slrdamage="mode")
    elseif AMOC=="T"
        m = base_model(rcp=Scenario=="ssp126" ? "RCP3-PD/2.6" : "RCP$(Scenario[5]).$(Scenario[6])", ssp="SSP$(Scenario[4])", tdamage=Dam=="BHM" ? "pointestimate" : Dam, slrdamage="mode")
        m = add_AMOC(m, Calib)
    elseif AMOC=="C"
        @unpack Eta, DecreaseCalib = d
        m = get_amoc_carbon_model(DecreaseCalib; scenario=Scenario, Dam=Dam)
        update_param!(m, :AMOC_Carbon_eta_AMOC, Eta)
    elseif AMOC=="T&C"
        println("Caution: this T&C version works with an untested temperature emulation and should not be used.")
        
        m = base_model(rcp="RCP4.5", ssp="SSP2", tdamage=Dam=="BHM" ? "pointestimate" : Dam, slrdamage="mode")
        replace!(m, :CO2Converter => CO2Converter_AMOC);
        Mimi.set_first_last!(m, :CO2Converter, first = 2010, last = 2200);
        set_param!(m, :CO2Converter, :co2_amoc, zeros(dim_count(model, :time)))
        run(m)
        
        m = add_AMOC_emulator(m)
        Mimi.set_first_last!(m, :AMOC_emulator, first = 2010, last = 2200);
        
        @unpack Eta, DecreaseCalib = d
        update_param!(m, :AMOC_emulator_eta_AMOC, Eta)
        
        if DecreaseCalib==false
            @unpack DecreaseAMOC = d
            update_param!(m, :AMOC_emulator_beta_AMOC, (DecreaseAMOC / (T_AT_2100 - T_AT_2010)).*ones(dim_count(m, :time)))
            update_param!(m, :AMOC_emulator, :T_AT_2010, T_AT_2010)
            update_param!(m, :AMOC_emulator, :AMOC_2010, 19.)
        else 
            m = calibrate_AMOC_decrease(m, DecreaseCalib)
        end
    else
        error("AMOC must be one of 'none', 'T', 'C', or 'T&C'")
    end
    update_param!(m, :Consumption, :damagepersist, persist);
    update_param!(m, :Utility_PRTP, prtp)
    run(m);
    return m
end

function run_amoc_mc_sccs(d::Dict; pulse_year=2020, pulse_size=10.)
    m = get_base_model_w_params(d)
    
    @unpack emuc, AMOC, MC_samples = d
    n_samples = length(MC_samples)
    
    if AMOC=="C" || AMOC=="T&C" 
        @unpack Eta, Eta_stderr, DecreaseCalib, Scenario = d
        if Scenario=="ssp245"
            std_AMOC_proj, std_AMOC_pi = get_proj_std_ssp245(DecreaseCalib)
        else
            std_AMOC_proj, std_AMOC_pi = 0., 0.
        end
        @time sccs = calculate_scc_base_mc(m, n_samples, false, false, false, pulse_year, pulse_size, emuc, calc_nationals=false, sample_id_subset=MC_samples,
        eta_AMOC_dist_scc=true, mean_eta_AMOC_scc=Eta, std_eta_AMOC_scc=Eta_stderr,
        AMOC_proj_dist_scc=true, std_AMOC_proj_scc=std_AMOC_proj,
        AMOC_pi_dist_scc=true, std_AMOC_pi_scc=std_AMOC_pi,
        )[:other]
    else
        @time sccs = calculate_scc_base_mc(m, n_samples, false, false, false, pulse_year, pulse_size, emuc, calc_nationals=false, sample_id_subset=MC_samples,
        eta_AMOC_dist_scc=false)[:other]
    end
    
    fulld = copy(d)
    fulld["sccs"] = sccs
    return fulld
end

function get_proj_std_ssp245(esm_calibration)
    esm_calibration = get_primary_model_name(esm_calibration, cmip6_model_names)
    amoc_ds = get_projection_ds(esm_calibration; scenario="ssp245")
    amoc_esm_std = amoc_ds["std"] |> Array
    mean_amoc_esm_std = mean(amoc_esm_std)
    
    AMOC_pi_vals = CSV.read("../data/CMIP6_amoc/AMOC_pi_values.csv", DataFrame)
    AMOC_pi_std = AMOC_pi_vals[:, esm_calibration][2]
    
    return mean_amoc_esm_std, AMOC_pi_std
end

function run_amoc_carbon_mc(d::Dict)
    @unpack AMOC, Dam, prtp, emuc, persist, Scenario, MC_samples = d
    n_samples = length(MC_samples)
    
    if AMOC=="none"
        m = base_model(rcp=Scenario=="ssp126" ? "RCP3-PD/2.6" : "RCP$(Scenario[5]).$(Scenario[6])", ssp="SSP$(Scenario[4])", tdamage=Dam=="BHM" ? "pointestimate" : Dam, slrdamage="mode")
        
        update_param!(m, :Consumption, :damagepersist, persist)
        update_param!(m, :Utility_EMUC, emuc)
        update_param!(m, :Utility_PRTP, prtp)
        run(m)
        
        @time mcres = sim_base(m, n_samples, false, false, false; 
        sample_id_subset=MC_samples)[:other]
        
    elseif AMOC=="C"
        @unpack Eta, Eta_stderr, DecreaseCalib = d
        m = get_amoc_carbon_model(DecreaseCalib; scenario=Scenario, Dam=Dam);
        update_param!(m, :AMOC_Carbon_eta_AMOC, Eta)
        
        update_param!(m, :Consumption, :damagepersist, persist)
        update_param!(m, :Utility_EMUC, emuc)
        update_param!(m, :Utility_PRTP, prtp)
        run(m)
        
        if Scenario=="ssp245"
            std_AMOC_proj, std_AMOC_pi = get_proj_std_ssp245(DecreaseCalib)
        else
            std_AMOC_proj, std_AMOC_pi = 0., 0.
        end
        
        @time mcres = sim_base(m, n_samples, false, false, false; 
        eta_AMOC_dist=true, mean_eta_AMOC=Eta, std_eta_AMOC=Eta_stderr,
        AMOC_proj_dist=true, std_AMOC_proj=std_AMOC_proj,
        AMOC_pi_dist=true, std_AMOC_pi=std_AMOC_pi,
        sample_id_subset=MC_samples,
        save_rvs=false,
        )[:other]
    else
        error("AMOC must be 'C' for this function")
    end
    
    fulld = copy(d)
    fulld["mcres"] = mcres
    return fulld

end