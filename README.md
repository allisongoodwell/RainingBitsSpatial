# RainingBitsSpatial
 codes and data files associated with Goodwell, JHM 2020 (submitted)

DOI of paper: submitted
DOI of dataset: TBD

This repository contains several python codes and results files for the paper
It's Raining Bits: patterns in directional precipitation persistence across the U.S.

Inputs to these codes are a folder of precipitation data in netcdf format, from the CPC gage-based daily gridded dataset, such that link to the folder should be updated as appropriate.  Some codes are set up to run in parallel, and number of pools should also be updated as appropriate depending on the computer used.  Some outputs of earlier files (thresholds, HXonly) will create input files for subsequent codes, so should be run in the order listed below:


Contents as follows:

Find_thresholds.py: detects minimum precipitation thresholds as described in main paper
Plot_thresholds.py: plots maps with seasonal thresholds
Find_HXonly.py: detects precipitation entropies only (was left out of main code)

RainingBits_annualanalysis_ITmeasures.py: one-year windows of IT measures
RainingBits_TrendsCorrs.py: uses results of one-year windows to detect trends and correlations

RainingBits_overallanalysis_ITmeasures.py: overall analysis of 70 year dataset
RainingBits_overallanalysis_makemaps.py: summarize overall analysis to obtain maps and statistics for all grid cells

