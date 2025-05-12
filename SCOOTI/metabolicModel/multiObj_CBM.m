function multiObj_CBM(config)
  addpath(genpath('utils'));
  % multiObj_CBM.m calls CFRinterface/DFAinterface to run CFR/DFA modeling and generate flux prediction based on transcriptomics/metabolomics data.
  %   This code is capable of predicting fluxes with GEMs Recon1/Recon2.2/Recon3D, multi-objectives, gene knockouts, and with/without pFBA.
  % Inputs:
  %   jj : integer, index number that is used to scan parameters. Default value is set to 1
  %   DFA_kappa : float, parameter to minimize flux activity coefficients
  %   CFR_kappa : float, parameter to reduce the upper bound of reactions associated with downregulated genes
  %   CFR_rho : float, parameter to increase the lower bound of reactions associated with upregulated genes
  %   COBRA_path : string, path to access the openCOBRA toolbox
  %   GEM_path : string, path to access the genome-scale metabolic models like Recon1
  %   model_name : string, name of the GEM, options including Recon1, Recon2.2, and Recon3D
  %   obj_candidate_list_file : string, path to access the candidate of multi-objectives, only applicable for single-objective models
  %   input_obj_tb : string, path to access the file of objective coefficients, only applicable for multiobjective models
  %   paraLen : integer, number to indicate how how many parameters that will be selected for scanning
  %   random_para : bool, enable random sampling of parameters in defined ranges
  %   init_objective : integer, the number is the index of objective list. In the default settings, 1 is no objective and 2 is biomass
  %   genekoflag : bool, enable in silico single-gene knock
  %   rxnkoflag : bool, enable in silico reaction knockout
  %   FSflag : bool, output flux sampling from solution space
  %   pfba : bool, minimizing sum of reaction fluxes or not
  %   medium_perturbation : bool, enable in silico single-metabolite depletion
  %   data_dir : string, path to access significant genes/proteins/metabolites as constraints
  %   prefix_name : string, name of experiment. it will be used for 
  %   medium : string, name of sources of extracellular metabolites. options including DMEMF12 and KSOM
  %   late_stage : string, name of late stage which also indicates the suffix of constraint files
  %   early_stage : string, name of early stage which also indicates the suffix of constraint files
  %   simulation : string, indicate which type of modeling method. options including CFR and DFA
  %   constraint : bool, enable adding constraints or not
  %   save_root_path : string, path to save the results of predicted fluxes
  %   CFR_model_path : string. path to access the constrained CFR models. This is only for stacking/multiomics integration
  %   pairwise_CFR_model : bool, match the name of the contrained CFR models or not
  %   extraWeight : float, apply to scale the importance of contrained CFR models for stacking different constraints
  %   algorithm : string, default value is iMAT but options including iMAT, MOOMIN, INIT (currently not working)
  %   data_series : string, default value is an empty string and the input paths should be separated by ","
  %   prefix_series : string, default value is an empty string and the input paths should be separated by ","
  %   medium_series : string, default value is an empty string and the input paths should be separated by ","
  % Output: 
  %   None; however, the predicted fluxes will be saved
 
%% Step 1: Set defaults and initialize
config = set_default_parameters(config);

if strcmp(config.prefix_name, 'model')
    config.simulation = 'CFR';
    config.dkappa = -1; config.ckappa = 1; config.crho = 1;
    config.constraint = 0; config.paraLen = 1; config.random_para = 0;
end

if ~iscell(config.data_series)
    [config.data_series, config.prefix_series, config.medium_series] = ...
        split_series(config.data_series, config.prefix_series, config.medium_series);
end

%% Step 2: Initialize COBRA
addpath(config.COBRA_path);
run initCobraToolbox(false);
changeCobraSolver('gurobi');

%% Step 3: Load GEM model
model_struct = load(config.GEM_path);
fn = fieldnames(model_struct);
config.model = model_struct.(fn{1});
if strcmp(config.model_name, 'Recon1')
    bio_obj_idx = find(contains(config.model.rxns, 'biomass_objective'));
    config.model.c(bio_obj_idx) = 1;
end

%% Step 4: Set objective candidates and weights
if isempty(config.obj_candidate_list_file)
    config.obj_candidate_list_file = './metabolicModel/GEMs/obj52_metabolites_shen2019.csv';
end
config.obj_candidates = readtable(config.obj_candidate_list_file);

if isempty(config.input_obj_tb)
    config.input_obj_tb = readtable('./metabolicModel/bufferFiles/buffer_objective_coefficients.csv');
    config.input_objective_weights = 0;
else
    config.input_obj_tb = readtable(config.input_obj_tb);
    config.input_objective_weights = 1;
end

%% Step 5: Constraint and data preprocessing
if config.constraint == 1
    if isempty(config.data_series) || all(cellfun(@isempty, config.data_series))
      disp('inininin')
        [config.data_series, config.prefix_series, config.medium_series] = ...
            batch_input_preprocess(config.data_dir, config.prefix_name, config.medium);
    end
    config.ctrl = 1;
else
    config.data_series = {'./metabolicModel/bufferFiles/buffer_constraint.xlsx'};
    config.prefix_series = {'2CBC_Israel'};
    config.medium_series = {config.medium};
    config.ctrl = 0;
end

%% Step 6: Load CFR models (if any)
CFR_models = load_CFR_models(config.CFR_model_path);

%% Step 7: Set simulation parameters (kappa/rho)
[config.kappa, config.rho] = determine_kappa_rho(config);

%% Step 8: Run simulation for each dataset
for data_idx = 1:length(config.data_series)
    config.data_idx = data_idx;
    if ~isempty(CFR_models)
      config.CFR_model = CFR_models{data_idx};
    else
      config.CFR_model = {};
    end

    if ~config.input_objective_weights & strcmp(config.prefix_name, 'model'),
      for ii=2:length(config.obj_candidates{:,2}),
        disp(ii)
        config.ii = ii;
        params = prepare_simulation_parameters(config);
        dispatch_simulation(config.simulation, params);
      end
    else
      config.ii = config.init_objective;
      params = prepare_simulation_parameters(config);
      dispatch_simulation(config.simulation, params);
    end
end
end % end for the function











