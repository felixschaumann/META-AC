using Test
include("../src/scc.jl")
include("../src/lib/AMOC_Carbon_functions.jl")

model_base = base_model(rcp="RCP4.5", ssp="SSP2", tdamage="pointestimate", slrdamage="mode")
update_param!(model_base, :Consumption, :damagepersist, 0.25)
run(model_base)
scc_base = calculate_scc(model_base, 2020, 10., 1.05)

model_base2 = full_model(rcp="RCP4.5", ssp="SSP2", saf=false, interaction=false, pcf=false, omh=false, amaz=false, gis=false, ais=false, ism=false, amoc=false)
update_param!(model_base2, :Consumption, :damagepersist, 0.25)
run(model_base2)
scc_base2 = calculate_scc(model_base2, 2020, 10., 1.05)

model_amoc_pattern = full_model(rcp="RCP4.5", ssp="SSP2", saf=false, interaction=false, pcf=false, omh=false, amaz=false, gis=false, ais=false, ism=false, amoc="IPSL")
p_vector = ones(dim_count(model, :time))
p_vector[266] = 0
update_param!(model_amoc_pattern, :AMOC_uniforms, p_vector)
update_param!(model_amoc_pattern, :Consumption, :damagepersist, 0.25)
run(model_amoc_pattern)
scc_amoc_pattern = calculate_scc(model_amoc_pattern, 2020, 10., 1.05)

model_amoc_nocarbon_pattern = get_amoc_carbon_model("mpiesmlr"; temp_pattern="IPSL")
update_param!(model_amoc_nocarbon_pattern, :AMOC_Carbon, :Delta_AMOC, 35)
update_param!(model_amoc_nocarbon_pattern, :AMOC_Carbon_eta_AMOC, 0)
run(model_amoc_nocarbon_pattern)
scc_amoc_nocarbon_pattern = calculate_scc(model_amoc_nocarbon_pattern, 2020, 10., 1.05)

model_amoc_nocarbon_nopattern = get_amoc_carbon_model("mpiesmlr"; temp_pattern=false)
update_param!(model_amoc_nocarbon_nopattern, :AMOC_Carbon_eta_AMOC, 0)
run(model_amoc_nocarbon_nopattern)
scc_amoc_nocarbon_nopattern = calculate_scc(model_amoc_nocarbon_nopattern, 2020, 10., 1.05)

@test scc_base == scc_base2
@test scc_base == scc_amoc_nocarbon_nopattern
@test scc_amoc_pattern == scc_amoc_nocarbon_pattern