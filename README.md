# mrb-streamflow-change
modeled changes in Missouri River flow due to climate, LULCC, and CO2 fertilization

## Script descriptions and order

### Estimate Q from model simulations

#### model_mstmip_contributions_runoff.m
Estimate runoff and streamflow for MsTMIP simulations 

#### model_climate_contributions_runoff.m
Estimate streamflow from water-budget model simulations driven by Williams variability vs. trends climate data

#### calculate_runoff_ratio.m
Compare modeled and measured runoff ratio (for supplemental material)

#### compare_models_gages.m
Compare modeled streamflow (MsTMIP and MW11) to measured streamflow at Hermann, MO gage (supplementary figure)

#### compare_wy_peak_mean.m
Compare water-year peak flow to water-year mean flow (supplementary figure)


### Make study area maps

#### map_study_area.m
Map of study area with major streams, regions, and instrumental and modeled flow anomalies

#### map_nlcd.m
Map of study area land cover (supplementary figure)

#### map_mstmip_lulcc.m
Map and time series of 1800-2010 land cover change using on the LULCC data used in the MsTMIP simulations (supplementary figure)


### Calculate driver contributions and make maps / tables of results

#### map_preindustrial_contributions_runoff.m
Calculate MsTMIP drivers of flow and make maps of contributions relative to "pre-Industrial" baselines

#### map_instrumental_contributions_runoff.m
Calculate drivers of flow (including anthropogenic vs. natural climate effects and separation of precipitaiton, PET, and temperature effects) and make maps of contributions to early instrumental baseline (i.e., earliest 30-year period of flow records)

#### make_table_dQ_PreIndustrial.m
Calculate contribution of each MsTMIP factor (climate, LULCC, and CO2) to flow change relative to "pre-Industrial" baselines

#### make_table_dQ_GageOverlap.m
Calculate contribution of each factor to flow change between mid and late century and compare to measured flow change

#### plot_monthly_climate_contributions_runoff.m
Calculate contributions of each factor to flow change between mid and late century for each individual region within the MRB

#### plot_qnat_teleconnection_correlations.m
Calculate and plot correlations between the "natural" component of MRB flow and seasonal composites of teleconnection indices (Nino3.4, PNA, NAO, and AMO)

