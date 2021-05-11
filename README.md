# Gradient Boosted Trees Capture Non-Linear Genetics Effects and Allele Interactions in Complex Phenotypes
## A Machine Learning Approach to the Polygenic Risk Score

This repository contains supporting code and documentation for the article "Gradient Boosted Trees Capture Non-Linear Genetics Effects and Allele Interactions in Complex Phenotypes".

The scripts in the pipeline must be run in the order specified, from a bash script. Each script builds on the last. This pipeline was run on a distributed cluster across 10 cores with 4 threads each.

## Documentation

`run_xgb_lasso_jointcv_v2.py`

Purpose: Cross-validate for optimal value of alpha. For each list of SNPs from the prior script, execute LASSO and XGBoost using only those SNPs. 

Inputs: Genotypes, phenotypes, non-zero coefficients from LASSO script 

Arguments:  

	phenotype: coded name of phenotype aligning to file names (options: TC, TG, SBP, SleepDur, Height) 

	var_phenotype: coded name of phenotype in the phenotype file (options: total_cholesterol_1, trigylcerides_1, bp_systolic_1, sleep_duration_1, height_baseline_1) 

	clumped_unclumped: whether to use the clumped or unclumped SNPs (options: unclumped or clumped) 

	rel: whether to include relatives in the training and validation set (options: Yes or No) 

	include_prs: whether to include PRS score in the algorithm (options: YesPRS or NoPRS) 

	race: race/ethnicity on which to train and assess algorithm (options: ALL, White, Black, HL) 

Outputs: For each alpha value, the validation RMSE and percentage of variance explained. You can use CV_Results_2021_0507.ipynb to identify the optimal values of alpha for each permutation of the above options. 

Output file names have the following format: 

'{}_{}_{}_{}Rel_{}_cvresults_{}.csv'.format(phenotype, race, clumped_unclumped, rel, include_prs, jobindex) 

------------------------------------------------------------------------------------------------------------------------ 

`run_xgb_lasso_v3_alpha.py`

Purpose: Given the optimal value of LASSO identified through cross-validation in the above script, execute and save the XGBoost models. 

Inputs: Genotypes, phenotypes, non-zero coefficients from LASSO script, optimal value of alpha 

Arguments:  

	phenotype: coded name of phenotype aligning to file names (options: TC, TG, SBP, SleepDur, Height) 

	var_phenotype: coded name of phenotype in the phenotype file (options: total_cholesterol_1, trigylcerides_1, bp_systolic_1, sleep_duration_1, height_baseline_1) 

	clumped_unclumped: whether to use the clumped or unclumped SNPs (options: unclumped or clumped) 

	rel: whether to include relatives in the training and validation set (options: Yes or No) 

	include_prs: whether to include PRS score in the algorithm (options: YesPRS or NoPRS) 

	race: race/ethnicity on which to train and assess algorithm (options: ALL, White, Black, HL) 

	alpha: optimal value of alpha with respect to the XGBoost objective function, identified from the prior analysis. 

Outputs: Model weights for the XGBoost model 

Output File names have the following format: 

'xgb_{}_{}_{}.sav'.format(var, race, clumped_unclumped) 

------------------------------------------------------------------------------------------------------------------------ 

`assess_xgb_newpop_v3.py`

Purpose: Given the saved XGBoost models for each ethnicity and optimal alpha values (identified in the 2 prior scripts), assess the XGBoost models on each permutation of ethnicity (e.g., trained on White, assessed on Black). 

Inputs: Genotypes, phenotypes, non-zero coefficients from LASSO script, optimal value of alpha, saved XGB models 

Arguments:  

	phenotype: coded name of phenotype aligning to file names (options: TC, TG, SBP, SleepDur, Height) 

	var_phenotype: coded name of phenotype in the phenotype file (options: total_cholesterol_1, trigylcerides_1, bp_systolic_1, sleep_duration_1, height_baseline_1) 

	clumped_unclumped: whether to use the clumped or unclumped SNPs (options: unclumped or clumped) 

	rel: whether to include relatives in the training and validation set (options: Yes or No) 

	include_prs: whether to include PRS score in the algorithm (options: YesPRS or NoPRS) 

	race: race/ethnicity on which to assess algorithm (options: ALL, White, Black, HL) 

	alpha: optimal value of alpha with respect to the XGBoost objective function, identified from the prior analysis. 

	train_race: the race the model was trained on (options: ALL, White, Black, HL) 

Outputs: Percentage of variance explained for the new population based on the models trained  using the trained population. 

Output File names have the following format: 

'assess_results_by_race_v3.txt' 


## Reference

Please reference this work as follows:



