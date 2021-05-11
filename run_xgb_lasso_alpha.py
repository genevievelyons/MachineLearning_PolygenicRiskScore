###############################################################
# Purpose :  Run XGB and LASSO for a given phenotype
# Created by : Genevieve Lyons
# Created date : 1/23/2021
# Update: Run by race and include PRS results 2/20/2021
###############################################################

import pandas as pd
import numpy as np
import glob
import os
import argparse

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression, Lasso, LassoCV
from sklearn.metrics import r2_score, explained_variance_score, mean_squared_error

import xgboost as xgb
import pickle

#######################################################################################################
# Parameters for Execution
#######################################################################################################

# Jobindex
jobindex = os.environ['LSB_JOBINDEX']
n_jobs = 5

# Parse arguments for phenotype
parser = argparse.ArgumentParser(description='XGB for given phenotype')
parser.add_argument('phenotype', type=str, help='phenotype')
parser.add_argument('var_phenotype', type=str, help='phenotype variable name')
parser.add_argument('clumped_unclumped', type=str, help='clumped or unclumped')
parser.add_argument('rel', type=str, help='include related individuals in training: Yes or No')
parser.add_argument('include_prs', type=str, help='include PRS results: YesPRS or NoPRS')
parser.add_argument('race', type=str, help='race: White, Black, HL, ALL')
parser.add_argument('alpha', type=str, help='alpha: chosen thru CV')
args = parser.parse_args()

phenotype = args.phenotype         # Add this in as an argument to the cmd line
var_pheno = args.var_phenotype         # Add this in as an argument to the cmd line
clumped_unclumped = args.clumped_unclumped         # Add this in as an argument to the cmd line
rel = args.rel        
include_prs = args.include_prs        
race = args.race         # Add this in as an argument to the cmd line
alpha = args.alpha         # Add this in as an argument to the cmd line

print(jobindex)
print(phenotype)
print(var_pheno)
print(clumped_unclumped)
print(race)
print(rel)
print(include_prs)

# Directories
model_weights_dir = ''
lasso_dir = ''

#######################################################################################################
## Function to load genotype and phenotype data ##
#######################################################################################################

def load_data(var, var_pheno, race, rel, lasso_lim=True, a=0, clumped_unclumped='unclumped'):
    '''
    Function to input a phenotype and alpha level, return the genotype and phenotype data

    Parameters:
        var : phenotype name in saved file, e.g. "TG" for triglycerides
        var_pheno : phenotype name in phenotype file, e.g. "triglycerides_1" for triglycerides
        lasso_lim : Whether to limit SNPs to those selected through LASSO feature selection (default True)
        a : alpha value to use for LASSO limitation
        clumped_unclumped : whether to use the clumped data through PRSice or the unclumped SNP data
        race : race to train and assess on, i.e. ALL, White, Black, or HL
        rel : whether to include 
        
    Returns:
        genotypes_train : SNP data for train set 
        genotypes_test : SNP data for test set
        y_train : target phenotype data for train set
        y_test : target phenotype data for test set
        phenotype_data_train : other covariate data for the train set (e.g., Age, Race, Sex)
        phenotype_data_test : other covariate data for the train set (e.g., Age, Race, Sex)
    '''
    
    ####################
    ## Genotype data ##
    ####################

    genotypes_train = []
    genotypes_test = []

    if lasso_lim==True:
        # Identify LASSO files for the alpha value
        fname_start = '{}_{}_SNPs_p1e-4*TOPMedTrain_Agefixed_MedicationFixedNoTCTGfix_Dedup{}Rel_genotypes_encoded_LASSO_part'.format(var, clumped_unclumped, rel)
        files = glob.glob(lasso_dir + fname_start + '*{}*Lambda_{}.csv.gz'.format(race, a)) 

        # Pull lasso limited data
        lasso = []
        for f in files:
            lasso.append(pd.read_csv(f))

        if len(lasso) >= 1:
            lasso = pd.concat(lasso, axis=0, ignore_index=True)
        else:
            print("Failed job - {} {} - job {} with alpha {}".format(var, clumped_unclumped, jobindex, a))
            lasso = lasso[0]

        # Parse the LASSO non-zero field names
        lasso = lasso.lasso_Names.str.split(';')
        # Stack the multiple series
        lasso = lasso.apply(pd.Series).stack().reset_index(drop=True)
        # Drop duplicates
        lasso = np.array(lasso.drop_duplicates())
        # Exclude intercept, reformat
        lasso = lasso[lasso!='(Intercept)']

    for c in range(1, 23):

        for cohort in ['TOPMedTrain', 'TOPMedTest']:
            # Define file name
            if cohort=='TopMedTrain':
                r = rel 
            else: 
                r = 'No'
            genotype_file = lasso_dir + '{}_{}_SNPs_p1e-4_{}_Agefixed_MedicationFixedNoTCTGfix_Dedup{}Rel_genotypes_encoded_chr_{}.csv.gz'.format(var, clumped_unclumped, cohort, r, c)
            # Load genotype data
            genotype = pd.read_csv(genotype_file)
            # Index as Subject ID
            genotype.index=genotype['Unnamed: 0']
            genotype.drop(['Unnamed: 0'], axis=1, inplace=True)

            if lasso_lim==True:
                # Only include the LASSO fields
                genotype = genotype[genotype.columns[genotype.columns.isin(lasso)]]

            # Append to genotype dataframe
            if cohort=='TOPMedTrain':
                genotypes_train.append(genotype)
            else:
                genotypes_test.append(genotype)

    # Flatten out the filtered genotype file
    genotypes_train = pd.concat(genotypes_train, axis=1, ignore_index=False, sort=False)
    genotypes_test = pd.concat(genotypes_test, axis=1, ignore_index=False, sort=False)
    
    ####################
    ## Phenotype data ##
    ####################
    # Load the phenotype data - train
    pheno_file = 'phenotypes_new_LIMITED_{}Rel_train_{}.csv'.format(rel, race)
    phenotype_data_train = pd.read_csv(pheno_file, dtype = {'SUBJECT_ID':'str', 'sample_id':'str'})
    # Index phenotype data on observation ID
    phenotype_data_train.index = phenotype_data_train.sample_id

    # Load the phenotype data - test
    pheno_file = 'phenotypes_new_LIMITED_validate_{}.csv'.format(race)
    phenotype_data_test = pd.read_csv(pheno_file, dtype = {'SUBJECT_ID':'str', 'sample_id':'str'})
    # Index phenotype data on observation ID
    phenotype_data_test.index = phenotype_data_test.sample_id
    
    ##########################
    ## Combine Data ##
    ##########################
    # Re-order genotype data as the pheno data - Train
    genotypes_train = genotypes_train.reindex(phenotype_data_train.index)

    # Re-order genotype data as the pheno data - Test
    genotypes_test = genotypes_test.reindex(phenotype_data_test.index)

    # Subset the pheno data to the EV fields
    y_train = phenotype_data_train[var_pheno]
    y_test = phenotype_data_test[var_pheno]

    phenotype_data_train = pd.concat([pd.get_dummies(phenotype_data_train[['Sex', 'race', 'study']]),
                                    phenotype_data_train[['EV1', 'EV2', 'EV3', 'EV4', 'EV5', 'Age', var_pheno]]],
                                  axis=1)
    phenotype_data_test = pd.concat([pd.get_dummies(phenotype_data_test[['Sex', 'race', 'study']]),
                                    phenotype_data_test[['EV1', 'EV2', 'EV3', 'EV4', 'EV5', 'Age', var_pheno]]], 
                                  axis=1)

    # Fix missing columns
    phenotype_data_test = phenotype_data_test.reindex(columns=phenotype_data_train.columns).fillna(0)

    # Exclude any subjects with NA's - train
    sna = phenotype_data_train.isna().any(axis = 1) | genotypes_train.isna().any(axis = 1)
    print("Subject with NA's - Train: {}".format(genotypes_train[sna == False].shape[0]))
    genotypes_train = genotypes_train[sna == False]
    y_train = y_train[sna == False]
    phenotype_data_train = phenotype_data_train[sna == False].drop(var_pheno, axis=1)

    # Exclude any subjects with NA's - test
    sna = phenotype_data_test.isna().any(axis = 1) | genotypes_test.isna().any(axis = 1)
    print("Subject with NA's - Test: {}".format(genotypes_test[sna == False].shape[0]))
    genotypes_test = genotypes_test[sna == False]
    y_test = y_test[sna == False]
    phenotype_data_test = phenotype_data_test[sna == False].drop(var_pheno, axis=1)

    # Exclude outliers using the 1st and 99th percentile
    genotypes_test = genotypes_test[(y_test >= y_train.quantile(0.01)) & (y_test <= y_train.quantile(0.99))]
    phenotype_data_test = phenotype_data_test[(y_test >= y_train.quantile(0.01)) & (y_test <= y_train.quantile(0.99))]
    y_test = y_test[(y_test >= y_train.quantile(0.01)) & (y_test <= y_train.quantile(0.99))]

    genotypes_train = genotypes_train[(y_train >= y_train.quantile(0.01)) & (y_train <= y_train.quantile(0.99))]
    phenotype_data_train = phenotype_data_train[(y_train >= y_train.quantile(0.01)) & (y_train <= y_train.quantile(0.99))]
    y_train = y_train[(y_train >= y_train.quantile(0.01)) & (y_train <= y_train.quantile(0.99))]

    # Print metrics
    print("n_train = {}".format(genotypes_train.shape[0]))
    print("n_test = {}".format(genotypes_test.shape[0]))
    print("Number of SNPs included - {}".format(genotypes_train.shape[1]))
    
    return genotypes_train, genotypes_test, y_train, y_test, phenotype_data_train, phenotype_data_test

#######################################################################################################
## Define function to Pull PRS Results ##
#######################################################################################################

def load_prsice(var, phenotype_data_train, phenotype_data_test, 
    clumped_dir = ''):
    
    # Load PRSice
    fname = clumped_dir + var + '_MAF0.01_R_0.1_250kb.all_score'
    prs = pd.read_csv(fname, sep = '\s+')

    # Set index
    prs.index=prs.FID
    prs.drop(['FID','IID'], axis=1, inplace=True)

    # Train/Test split
    prs_train = prs.reindex(phenotype_data_train.index)
    prs_test = prs.reindex(phenotype_data_test.index)
    
    return prs_train, prs_test

#######################################################################################################
## Define function to run linear model on age and sex and return the adjusted phenotype ##
#######################################################################################################

def adj_pheno(pheno_train, pheno_test, y_train, y_test):
    
    pheno_train = pheno_train.drop(['EV1', 'EV2', 'EV3', 'EV4', 'EV5'], axis=1)
    pheno_test = pheno_test.drop(['EV1', 'EV2', 'EV3', 'EV4', 'EV5'], axis=1)
    
    # Run lm
    lm = LinearRegression(fit_intercept=True, n_jobs=-1).fit(pheno_train, y_train)

    # Calculate genetic component
    y_train_pred = lm.predict(pheno_train)
    y_test_pred = lm.predict(pheno_test)

    # Print Metrics
    print("Train R^2 Covariates: {}".format(r2_score(y_train, y_train_pred)))
    print("Test R^2 Covariates: {}".format(r2_score(y_test, y_test_pred)))
    
    print("Train EVR Covariates: {}".format(explained_variance_score(y_train, y_train_pred)))
    print("Test EVR Covariates: {}".format(explained_variance_score(y_test, y_test_pred)))
    
    # Calculate residuals
    y_train_resid = y_train - y_train_pred
    y_test_resid = y_test - y_test_pred
    
    return y_train_resid, y_test_resid

#######################################################################################################
## Define function to run LASSO, cross validate for the optimal number of trees, and report out metrics ##
#######################################################################################################

def run_lasso(var, X_train, X_test, y_train, y_test, save=True):

    # Run LASSO with CV
    lasso = LassoCV(cv = 3, fit_intercept = True, n_jobs = -1, normalize = False).fit(X_train, y_train)

    # Save model
    if save==True:
        pickle.dump(lasso, open(model_weights_dir + 'lasso_{}.sav'.format(var), 'wb'))
        # Example: load the model from disk
        # lasso = pickle.load(open(model_weights_dir + 'lasso_{}.sav'.format(var), 'rb'))

    # Calculate genetic component
    genetic_component_train = lasso.predict(X_train)
    genetic_component_test = lasso.predict(X_test)

    # Print Metrics
    # print("Train R^2 LASSO: {}".format(r2_score(y_train, genetic_component_train)))
    # print("Test R^2 LASSO: {}".format(r2_score(y_test, genetic_component_test)))
    
    print("Train EVR LASSO: {}".format(explained_variance_score(y_train, genetic_component_train)))
    print("Test EVR LASSO: {}".format(explained_variance_score(y_test, genetic_component_test)))

    test_evr = explained_variance_score(y_test, genetic_component_test)
    test_rmse = np.sqrt(mean_squared_error(y_test, genetic_component_test))
    
    return genetic_component_train, genetic_component_test, test_rmse, test_evr

#######################################################################################################
## Define function to run XGB, cross validate for the optimal number of trees, and report out metrics ##
#######################################################################################################

def run_xgb(var, X_train, X_test, y_train, y_test, 
            params={'obj': 'reg:squarederror',
                    'max_depth': 5,
                    'colsample_bytree': 0.9,
                    'eta': 0.01,
                    'alpha': 0,
                    'gamma': 0,
                    'min_child_weight': 10,
                    'subsample': 0.5,
                    'nthread': -1},
            save=True):

    # Convert to XGBoost D-Matrix
    D_train = xgb.DMatrix(X_train, label=list(y_train))
    D_test = xgb.DMatrix(X_test, label=list(y_test))

    # Run XGBoost with CV for number of trees
    cv_results = xgb.cv(params, 
                        D_train,  
                        num_boost_round=10000, 
                        early_stopping_rounds=10, 
                        nfold=3,
                        as_pandas=True,
                        metrics={'rmse'})
    
    best_B = cv_results['test-rmse-mean'].idxmin()
    mean_mae = cv_results['test-rmse-mean'].min()
    print("RMSE {} with Best B = {}".format(mean_mae, best_B))
    
    # Train model with best B from CV
    xgb_model = xgb.train(params, D_train, num_boost_round = best_B)

    # Save model
    if save==True:
        pickle.dump(xgb_model, open(model_weights_dir + 'xgb_{}_{}_{}.sav'.format(var, race, clumped_unclumped), 'wb'))
        # Example: load the model from disk
        # xgb_model = pickle.load(open(model_weights_dir + 'xgb_{}.sav'.format(var), 'rb'))

    # Calculate genetic component
    genetic_component_train = xgb_model.predict(D_train)
    genetic_component_test = xgb_model.predict(D_test)

    # Print Metrics
    # print("Train R^2 XGB: {}".format(r2_score(y_train, genetic_component_train)))
    # print("Test R^2 XGB: {}".format(r2_score(y_test, genetic_component_test)))
    
    print("Train EVR XGB: {}".format(explained_variance_score(y_train, genetic_component_train)))
    print("Test EVR XGB: {}".format(explained_variance_score(y_test, genetic_component_test)))
    
    test_evr = explained_variance_score(y_test, genetic_component_test)
    
    return genetic_component_train, genetic_component_test, mean_mae, test_evr

#######################################################################################################
#######################################################################################################

########################################
# Run for pre-chosen optimal Alpha
########################################

# Load Data
X_train, X_test, y_train, y_test, pheno_train, pheno_test = load_data(var=phenotype, 
            var_pheno=var_pheno, race=race, rel=rel, lasso_lim=True, a=alpha,
            clumped_unclumped=clumped_unclumped)

print("Data Load Complete")

# Load PRSice
prs_train, prs_test = load_prsice(var=phenotype, phenotype_data_train=pheno_train, phenotype_data_test=pheno_test)

# Adjust the phenotype for age and sex
y_train, y_test = adj_pheno(pheno_train, pheno_test, y_train, y_test)

if include_prs=='NoPRS':
    train_df = pd.concat([X_train, pheno_train], axis=1)
    test_df = pd.concat([X_test, pheno_test], axis=1)
else:
    train_df = pd.concat([X_train, pheno_train, prs_train], axis=1)
    test_df = pd.concat([X_test, pheno_test, prs_test], axis=1)

# Run XGB 
params={'obj': 'reg:squarederror',
                        'max_depth': 5, # 5
                        'colsample_bytree': 0.9,
                        'eta': 0.01,
                        'alpha': 0,
                        'gamma': 0,
                        'min_child_weight': 10,
                        'subsample': 0.9, # 0.5
                        'nthread': -1}

genetic_component_train, genetic_component_test, mean_rmse_cv_xgb, test_evr_xgb = run_xgb(phenotype, 
                train_df, test_df, 
                y_train, y_test, 
                params=params,
                save=True)


