# ThesisCode

The folder for the thesis project.

## data

This folder contains the preprocessed datasets in ".Rda" format. This repository already contains the preprocessed datasets when it is downloaded, so the preprocessing does not have to be done in order to reproduce the results. The preprocessing can be done manually if desired by running the scripts in the "preprocessing" directory. 
 
## preprocessing

The folder for preprocessing codes which generates results in the folder "data/". The original datasets for the metz and davis dataset will automatically be downloaded from "http://staff.cs.utu.fi/~aatapa/data/DrugTarget/known_drug-target_interaction_affinities_pKi__Metz_et_al.2011.txt" and "http://staff.cs.utu.fi/~aatapa/data/DrugTarget/drug-target_interaction_affinities_Kd__Davis_et_al.2011.txt", when running the script "preprocess_metz_davis.R". The Kiba dataset must be downloaded manually from "http://pubs.acs.org/doi/suppl/10.1021/ci400709d". The filename that needs to be downloaded in order to run "preprocess_Kiba.R", is "ci400709d_si_002.xlsx". Note that it is not necessary to run the preprocessing step in order to reproduce the results, because the "data/" folder already contains the preprocessed datasets. 

## mf

The folder for the matrix factorization method

## kernel

The folder for the kernel method. This folder contains the code for the comparison model KronRLS. 

## xgboost

The folder for the xgboost method. This folder contains the code for the SimBoost and SimBoostQuant methods. See the README file in this folder for further instructions on how to use the method. 

## evaluation

The folder for the evaluation code. See the README file in this folder for further instructions.

