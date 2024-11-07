SCOOTI Settings
===============

Overview
--------
SCOOTI includes three different modules, including identifying metabolic regulators, flux predictions, and inference of metabolic objectives. This document will go through the required settings for users to choose the best options for their projects.

Identification of metabolic regulators
--------------------------------------

.. code-block:: bash

   #!/bin/bash


Flux prediction
---------------

**Usage:**

.. code-block:: bash

   #!/bin/bash

   matlab -nosplash -noFigureWindows -r "multiObj_CBM($run, $DFA_kappa, $CFR_kappa, $CFR_rho, $COBRA_path,$GEM_path, $model_name, $obj_candidate_list_file, $input_obj_tb, $paraLen, $random_para, $init_objective, $genekoflag, $rxnkoflag, $FVAflag, $pfba, $medium_perturbation, $data_dir, $prefix_name, $medium, $late_stage, $early_stage, $simulation, $constraint, $save_root_path)"

Input settings
**************

*multiObj_CBM.m calls CFRinterface/DFAinterface to run CFR/DFA modeling and generate flux prediction based on transcriptomics/metabolomics data.*
*This code is capable of predicting fluxes with GEMs Recon1/Recon2.2/Recon3D, multi-objectives, gene knockouts, and with/without pFBA.*

**Inputs:**
  - *jj* : integer, index number that is used to scan parameters
  - *DFA_kappa* : float, parameter to minimize flux activity coefficients
  - *CFR_kappa* : float, parameter to reduce the upper bound of reactions associated with downregulated genes
  - *CFR_rho* : float, parameter to increase the lower bound of reactions associated with upregulated genes
  - *COBRA_path* : string, path to access the openCOBRA toolbox
  - *GEM_path* : string, path to access the genome-scale metabolic models like Recon1
  - *model_name* : string, name of the GEM, options including Recon1, Recon2.2, and Recon3D
  - *obj_candidate_list_file* : string, path to access the candidate of multi-objectives, only applicable for single-objective models
  - *input_obj_tb* : string, path to access the file of objective coefficients, only applicable for multiobjective models
  - *paraLen* : integer, number to indicate how how many parameters that will be selected for scanning
  - *random_para* : bool, enable random sampling of parameters in defined ranges
  - *init_objective* : integer, the number is the index of objective list. In the default settings, 1 is no objective and 2 is biomass
  - *genekoflag* : bool, enable in silico single-gene knock
  - *rxnkoflag* : bool, enable in silico reaction knockout
  - *FVAflag* : bool, output a range of flux solution
  - *pfba* : bool, minimizing sum of reaction fluxes or not
  - *medium_perturbation* : bool, enable in silico single-metabolite depletion
  - *data_dir* : string, path to access significant genes/proteins/metabolites as constraints
  - *prefix_name* : string, name of experiment. it will be used for 
  - *medium* : string, name of sources of extracellular metabolites. options including DMEMF12 and KSOM
  - *late_stage* : string, name of late stage which also indicates the suffix of constraint files
  - *early_stage* : string, name of early stage which also indicates the suffix of constraint files
  - *simulation* : string, indicate which type of modeling method. options including CFR and DFA
  - *constraint* : bool, enable adding constraints or not
  - *save_root_path* : string, path to save the results of predicted fluxes
  - *CFR_model_path* : string. path to access the constrained CFR models. This is only for stacking/multiomics integration
  - *pairwise_CFR_model* : bool, match the name of the contrained CFR models or not
  - *extraWeight* : float, apply to scale the importance of contrained CFR models for stacking different constraints
  - *algorithm* : string, default value is iMAT but options including iMAT, MOOMIN, INIT (currently not working)
  - *data_series* : string, default value is an empty string and the input paths should be separated by ","
  - *prefix_series* : string, default value is an empty string and the input paths should be separated by ","
  - *medium_series* : string, default value is an empty string and the input paths should be separated by ","

**Output:**
  - *None*; however, the predicted fluxes will be saved


Objective Inference
-------------------

**Usage:**

.. code-block:: python

   #!/usr/bin/env python3

   python SCOOTI_trainer.py [-h] [--unconModel] [--conModel] [--savePath] [--kappaArr] [--rhoArr] [--dkappaArr] [--expName] [--unconNorm] [--conNorm] [--medium] [--method] [--model] [--inputType] [--clusterPath] [--objListPath] [--rank] [--stackModel] [--sampling] [--learner] [--geneKO] [--geneListPath] [--learningRate] [--epo]

Input settings
**************

*SCOOTI_trainer.py calls regressoirTraining.py to run meta-learner regression with unconstrained models as variables and constrained models as outcomes to infer metabolic objectives.*
*This code is capable of inferring metabolic objectives with linear regression or ADALINE neuron as the meta-learner.*

**Inputs settings**

  -**--unconModel** : path to access the directory that saves unconstrained models (optimal models)
  -**--conModel** : path to access the directory that saves the constrained models (context-specific models)
  -**--savePath** : path to save the output of metabolic objectives
  -**--kappaArr** : parameter array of kappa for CFR
  -**--rhoArr** : parameter array of rho for CFR
  -**--dkappaArr** : parameter array of kappa for DFA
  -**--expName** : name of experiment used to save the objectives
  -**--unconNorm** : normalization of unconstrained models
  -**--conNorm** : normalization of constrained models
  -**--medium** : concentration of extracellular metabolites
  -**--method** : options including `cfr`, `dfa`, and `model`
  -**--model** : genome-scale metabolic model. Options including `recon1`, `recon2`, and `recon3`
  -**--inputType** : different types of inputs. Options including `flux` and `pageRank`
  -**--clusterPath** : path to access the csv file that saves the clusters of metabolites
  -**--objListPath** : path to access the csv file that saves the list of objective candidates
  -**--rank** : T or F, fit the meta-learner models with ranks of fluxes (input values)
  -**--stackModel** : T or F, reuse the CFR models and constrain the models with new omics data
  -**--sampling** : T or F, fit the meta-learner models with sampled fluxes
  -**--learner** : A or L, apply either linear regressor or ADALINE as the meta-learner model
  -**--geneKO** : T or F, enable modeling single gene knockouts
  -**--geneListPath** : path to access the csv file that save the list of genes
  -**--learningRate** : learning rate aka step size for ADALINE neuron
  -**--epo** : limit of learning cycles for ADALINE neuron



