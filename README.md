# mrb-streamflow-change
modeled changes in Missouri River flow due to climate, LULCC, and CO2 fertilization

## Script descriptions and order

### Estimate Q from model simulations

#### model_mstmip_contributions_runoff.m
Estimate runoff and streamflow for MsTMIP simulations (and rescale to match mean and variance of gage streamflow)

#### model_climate_contributions_runoff.m
Estimate streamflow from water-budget model simulations driven by Williams variability vs. trends climate data

#### calculate_runoff_ratio.m
Compare modeled and measured runoff ratio (for supplemental material)

#### compare_models_gages.m
Compared modeled streamflow (unscaled MsTMIP, scaled MsTMIP, and WBM) to measured streamflow at Hermann, MO gage


### Make study area maps

#### map_study_area.m
Map of study area with major streams, regions, and instrumental flow

#### map_nlcd.m
Map of study area land cover (for supplemental material)

#### map_mstmip_lulcc.m
Map and time series of 1800-2010 land cover change using on the LULCC data used in the MsTMIP simulations


### Calculate driver contributions and make maps / tables of results

#### map_mstmip_contributions_runoff.m
Calculate MsTMIP drivers of flow and make maps of contributions

#### map_climate_contributions_runoff.m
Calculate interannual variability and anthropogenic trend contributions to flow (including separation of precipitaiton, PET, and temperature) and make maps of contributions

#### make_table_dQ.m
Calculate contribution of each factor to flow change between mid and late century and compare to measured flow change

#### make_mstmip_contributions_table.m
Make table of contributions from MsTMIP drivers (relative to MsTMIP baselines) to flow change for each subregion

