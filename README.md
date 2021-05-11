# Steady State Evoked Potential (SSVEP)

This repository contains code for the CRAMPP Visual Task EEG processing and analysis.

# Order of Processing

The following describes the order of processing the EEG files (i.e., pipeline) and important features of this workflow.

Two important files/scripts that nearly all EEG scripts (`*.m`) require are:

* `src/vis_workspace_prep.m` -> prepares workspace by setting working directory, starting EEGlab, and reading in participant information

* `data/0-global-vars/vis-subj-info.xlsx'` -> this workbook stores notes about the subjects (e.g., data quality, why participants were dropped, etc.). Most importantly, the first sheet serves as a batch processor in which one or more participants can be added for batch processing

Useful scripts/tools for preparing data:

* `rename_brainvision_files.m` -> script written by Robert Oostenveld (click here for [gist](https://gist.github.com/robertoostenveld/e31637a777c514bf1e86272e1092316e)) that helps rename brainvision files. If a brainvision EEG file is re-named after the recording is finished, then it becomes unreadable because the header expects the original name.

Preprocessing pipeline:

1. `prepro.m` -> preprocesses EEG files in this order:

* Channel locations are added
* Averged mastoid reference computed
* Downsamples data to 256Hz
* Removes DC offset
* Removes line noise using Cleanline plugin

2. Data are then visually inspected and bad channels (e.g., channels that were disconnected during recording) are identified and noted in `data/0-global-vars/vis-subj-info.xlsx'`

3. `ica_mara.m` -> performs Independent Components Analysis (ICA) and identifies artifactual ICs using [MARA](https://github.com/irenne/MARA)

ICA is performed on the entire recording to provide MARA a list of ICs to classify as artifactual. Importantly, ICA is not performed on bad channels and these are removed for MARA. Also, MARA's threshold was adjusted to be more conservative from .5 {default} to .6. All these results are saved out for each individual participant.

4. `epoching_mara.m` -> cleans and epochs data

* Artifactual ICs are dropped
* Bad channels are interpolated
* Stimulation blocks are extracted: -7 seconds before stimulation starts (to acquire a baseline) to 22 seconds after stimulation begins. 
* Within these new epochs, smaller 2 second epochs are created, with 1 second overlap (this is to overcome stationarity violations in FFTs)
* Epochs are then thresholded +- 100 microvolts and bad epochs are rejected

5. `spec_power_mara.m` -> spectral power calculations with surface Laplacian filter

Data are first spatially filtered using a surface laplacian (CSDtoolbox) and then filtered data are subject to `spectopo()` that computes spectral power (dB) at .5Hz bins

# Data Analysis and Visualization Scripts

Data analysis and visualization was performed in R.

1. Preprocessing scripts

* `eprime-prepro.R` -> Prepares visual task E-Prime (behavioral) data
* `medical-data-prepro.R` -> Prepares medical data (demographics etc.) 
* `spectral-prepro-dB.R`-> Prepares final spectral results from MATLAB/EEGLAB
* `ss-info-prepro.R` -> Prepares subject info (group classification etc.)


2. Useful scripts
* `topo_tools.R` -> Helpful functions for visualizing EEG topographies

3. Analysis scripts
* `analysis-8d.R`-> main analysis (multilevel modeling, visualization, etc.)
* `analysis-8d-supplement.R`-> supplementary analysis
