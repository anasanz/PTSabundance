# Data and code for manuscript “How to better count elusive birds? Comparing non-invasive monitoring methods to estimate population size of the endangered Pin-tailed sandgrouse (Pterocles alchata)” 

## Description of the data and file structure

### Code
This repository contains the following code files:

#### <ins>GMR folder</ins>
- *CR_model_run_results.r*: R code to arrange the data, run the GMR model and process results
- *simulation_test_accuracy.r*: R code of the simulation performed to test the accuracy of the GMR approach under our study conditions  

<ins>optimize_cost folder</ins>: analyses to optimize costs of GMR approach. This folder includes:
  - *1.optimize_subsample_dataSum*: Run CMR model by randomly selecting between 50% and 95% of the samples with successful genetic identification. Results are stored in *results_per2.RData*, which are provided as data files and loaded in script 3   
  - *2.optimize_remove_oc*: Run CMR model by removing between 1 and 3 sampling occasions in all possible combinations. Results are stored in *results_oc_remove.RData*, which are provided as data files and loaded in script 3  
  - *3. process_results*: Assess bias and CV from cost optimization analyses (scripts 1 and 2)  

#### <ins>GeneralHDS</ins>
- *MultispeciesHDS_model_run_results.r*: R code to arrange the data, run the MultispeciesHDS model and process results. Model results are stored in *3.Dat_HDS_trendmodel_lam[hq_weight_CAT2]_sigHR.RData*, which are provided as data files  
- *3.HDS_trendmodel_lam[hqCAT]_sigHR_noW*: MultispeciesHDS model  
  
#### <ins>SpecificHDS</ins>
- *SandgrouseHDS_model_run_results.r*: R code to arrange the data, run the SandgrouseHDS model and process results.
- *1.2.HDS_sig[HR_inf_fullModel]_lam[hq[inf]_CAT]_corrected*: SpecificHDS model
- *simulation_repeatedSurveys_sandgrouseHDS*: R code for simulation performed to test how many HDS surveys would be needed to reach the precision of GMR
- *5.HDS_2Survey_sig[HR_inf_fullModel]_lam[hq[inf]_CAT]*: SpecificHDS model with additional survey dimension for the simulation performed in the script *simulation_repeatedSurveys_sandgrouseHDS*

### Data
The Zenodo repository [https://zenodo.org/records/14387266?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6Ijc5NjlkYTdiLTY0ZTItNDVjZi05OGU1LWFiYzFjZTc3ODQxYSIsImRhdGEiOnt9LCJyYW5kb20iOiI3MDRiMzllYWYyMjJlNjU1ZjBjMWUxMDk2YzhlYjlkMSJ9.-XL0ckizaV1vrU6spF7Pjh0iVvmluqFyM3Hz_mDNLu7LssAymPSepS-QKvdYYBGxFhFH5pzO7f8jxG1IIP9Vdw](https://zenodo.org/records/14387266?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6Ijc5NjlkYTdiLTY0ZTItNDVjZi05OGU1LWFiYzFjZTc3ODQxYSIsImRhdGEiOnt9LCJyYW5kb20iOiI3MDRiMzllYWYyMjJlNjU1ZjBjMWUxMDk2YzhlYjlkMSJ9.-XL0ckizaV1vrU6spF7Pjh0iVvmluqFyM3Hz_mDNLu7LssAymPSepS-QKvdYYBGxFhFH5pzO7f8jxG1IIP9Vdw) contains all data files needed to run the aboved described scripts. Please ask corresponding author to get access to the data files.
