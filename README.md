# simus_orca1_fluxforced

Force the ocean in NEMO3.6 (no sea ice) with fixed fluxes at the ocean surface. Configuration based on the ORCA1_LIM3_PISCES configuration. Here you will find the modified code to put in the MY_SRC directory, as well as the namelist_cfg for the forced simulations. Detailed instructions concerning the modifications in each routine can be found in the document "Forcage_en_flux_ORCA1" (in French). The implementation of two passive tracers is also included (trc_sms_mytrc).

The configuration was developped in the context of my PhD work, with the help of Clément Rousset, Christian Ethé and Gurvan Madec (LOCEAN-IPSL). 
