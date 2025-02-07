### SRM in a nutshell

Code to post-process common SAI simulations into a small (~3GB), user-friendly, archive of data in an analysis-ready format.


Inputs:
* ARISE-1.5K simulations (UKESM1 and CESM2-WACCM)
* G6sulfur simulations for 6 models: UKESM1, CESM2-WACCM, IPSL-CM6A, MPI-ESM1-2 (High and low resolution), CESM2-WACCM, and CNRM-ESM2

For each SAI simulation, background/control runs are also used. For ARISE, this is the SSP2-4.5 run (and the end of the historical run for UKESM1). For G6sulfur, this is both the SSP5-8.5 run (the 'background') and the SSP2-4.5 run (the 'baseline' and 'target').

Outputs:



To run the code in Usage_example/, first download the data archive from zenodo (https://zenodo.org/records/14802397), then place this unzipped archive in a folder 'pp_archive/' at the root of this project directory.  



A note on calendars: 
UKESM uses a 360 day calendar, so there is no need to apply month weights when calculating annual means. For CESM2 and others, we do need to do this, hence the additional step in meaning in these scripts, using month weighting functions (in utils.py).


