function multiObj_CBM(jj, DFA_kappa, CFR_kappa, CFR_rho, COBRA_path, GEM_path, model_name, obj_candidate_list_file, input_obj_tb, paraLen, random_para, init_objective, genekoflag, rxnkoflag, FVAflag, pfba, medium_perturbation, data_dir, prefix_name, medium, late_stage, early_stage, simulation, constraint, save_root_path)
  % multiObj_CBM.m calls CFRinterface/DFAinterface to run CFR/DFA modeling and generate flux prediction based on transcriptomics/metabolomics data.
  %   This code is capable of predicting fluxes with GEMs Recon1/Recon2.2/Recon3D, multi-objectives, gene knockouts, and with/without pFBA.
  % Inputs:
  %   jj : integer, index number that is used to scan parameters
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
  %   FVAflag : bool, output a range of flux solution
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
  % Output: 
  %   None; however, the predicted fluxes will be saved
 
%% default settings for DFA and CFR parameters
  if (~exist('DFA_kappa','var')) || (isempty(DFA_kappa))
      dkappa = 1;
  end
  
  if (~exist('CFR_kappa','var')) || (isempty(CFR_kappa))
      ckappa = 1;
  end
  if (~exist('CFR_rho','var')) || (isempty(CFR_rho))
      crho = 1;
  end
%% Part 1: metabolic model
  % initialize COBRA toolbox
  addpath(COBRA_path)
  run initCobraToolbox(false)
  % INPUT files
  % Load metabolic network model
  model = load(model_path);
  fn = fieldnames(model);
  model = getfield(model, fn{1});
  % fix the model for the default objective
  if strcmp(model_name, 'Recon1'),
    model.c(find(contains(model.rxns, 'biomass_objective'))) = 1;
  end
  % objective candidates
  obj_candidate_list_file = '';

% Part 2: Input data
  % random sampling
  reind = randsample(paraLen*paraLen, paraLen);
  input_data_choose = 'model';

  % additional input for objectives and corresponding weights
  if strcmp(input_obj_tb, ''),
    input_obj_file = './metabolicModel/buffer_objective_coefficients.csv';
    input_objective_weights = 0; % set to 1 if giving existing objective weights
  else,
    input_obj_file = input_obj_tb;
    input_objective_weights = 1;
  end
  
  % read coefficients for objectives
  input_obj_tb = readtable(input_obj_file);

  %% Significant genes/proteins/metabolites
  if constraint==1,
    % data input setting
    [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);
    ctrl = 1 %set to 1 if we want to apply constraints
  else, % CFR without constraint
    data_series = {'./metabolicModel/buffer_constraint.xlsx'};
    prefix_name = 'unconstraint';
    prefix_series = {'2CBC_Israel'};
    ctrl = 0; % remove constraint
    medium_series = {medium};
  end


%% Part 3: Parameter settings for CBM 
  if strcmp(simulation, 'DFA')==1, % optimize reactions/new demand reactions
%% settings for parameter space for DFA kappa
    if dkappa==-1,
      DFA_para_space = logspace(1, -3, paraLen);
      disp('Parameter settings for DFA...')
      if random_para==1,
        DFA_para_space = logspace(1, -3, paraLen*paraLen);
        DFA_para_space = DFA_para_space(reind);
      end
      kappa = DFA_para_space(jj)
      rho = 0; % for bugs in DFA (unnecessary parameter)
    else,
      kappa = dkappa;
      rho = 0;
    end
  else, % CFR or model
%% settings for parameter space of CFR kappa and rho
    if ckappa==-1,
      CFR_para_space = logspace(1, -3, paraLen);
      [Y, Z] = meshgrid(CFR_para_space, CFR_para_space);
      CFR_para_pack = [Y(:), Z(:)];
      disp('Parameter settings for CFR...')
      if random_para==1,
        CFR_para_pack = CFR_para_pack(reind, :);
      end
      kappa = CFR_para_pack(jj, 1);
      rho = CFR_para_pack(jj, 2);
    else,
      kappa = ckappa;
      rho = crho;
    end
  end

%% Part4: Iterate through all different input datasets
  % setups for parallel computing
  delete(gcp('nocreate'));
  parpool(32);
  %disp(length(data_series))
%% Parfor...
  parfor data_index=1:length(data_series),
  %for data_index=1:length(data_series),
    %data_index = 1;
    size(data_series)
    size(medium_series)
    data_path = data_series{data_index}
    disp('Read file...')
    disp(data_path)
    % info of the data/simulation
    prefix = prefix_series{data_index}
    % get cell-specific objective with corresponding cell names
    if input_objective_weights & strcmp(input_data_choose, 'model')==0,
      prefix_str = prefix_pattern_func(prefix);
      input_obj_tb_ind = find([logical(1) contains(input_obj_tb.Properties.VariableNames(2:end), prefix_str)]);
    end
    % change culture medium KSOM_AA_Jin or DMEMF12_MGSA
    medium = medium_series{data_index};
    
    % candidate objectives
    single_obj = {};
    single_obj{1, 1} = '';
    candidates = readtable(obj_candidate_list_file);
    candidate_size = size(candidates);
    for i=1:candidate_size,
      single_obj(i+1, 1) = candidates{i, 2};
    end
    
    % objective weights for single objectives
    para_pack = eye(length(single_obj));
   
%% CFR simulations with different objective weights
    % parallel computing settings
    %init_for = 1;
    disp('input_obj_check')
    disp(input_objective_weights)
    if input_objective_weights==0, % single objective
      init_for = init_objective;
      if strcmp(input_data_choose, 'model')==1, % iterate thru all obj
        end_for = length(single_obj);
      else, % without obj for constraint models
        end_for = init_for;%length(single_obj);%1
      end
    elseif strcmp(input_data_choose, 'model')==1,
      % use new objectives without constraint
      init_for = 2;
      ss = size(input_obj_tb)
      end_for = ss(2);
    else,
      % use new objective WITH constraint
      init_for = 2;%init_for;
      ss = size(input_obj_tb(1:end, input_obj_tb_ind))
      end_for = ss(2);
    end
    %delete(gcp('nocreate'))
    %parpool(32);
    for ii=init_for:end_for,
    %parfor ii=init_for:end_for,

      % use single objective functions
      if input_objective_weights==0,
        % for single objective
        obj = single_obj; % comment out this one if you to model
        % objective weights
        disp('Weights for objectives:')
        obj_c = [];
        for d=1:length(obj),
          obj_c = [obj_c, para_pack(ii, d)];
        end

      % use multiobjective functions from files without constraint
      elseif strcmp(input_data_choose, 'model')==1,
        
        obj = input_obj_tb{2:ss(1), 1};
        obj_c = input_obj_tb{2:ss(1), ii}; 
        get_late_stage = input_obj_tb.Properties.VariableNames{ii};

      % use multiobjective function from files WITH corresponding constraint
      else,
        obj = input_obj_tb{2:ss(1), 1}
        obj_c = input_obj_tb{2:ss(1), ii}
        get_late_stage = late_stage

      end
      % DO NOT CHANGE ANYTHING BELOW
      % Running exps for all conditions
      disp('Optimizing the demand reaction:')
      disp(obj)
      if strcmp(obj, ''),
          obj_type = '';
      else,
          obj_type = 'Demand';
      end
      
      
      %%%%%%%%%%%%%%%%%%
      % start simulation
      %%%%%%%%%%%%%%%%%%
      % later stages
      
      % filename logics: 
      % 1: '[datetime]Project_simulation_fluxes.csv'
      % 2: '[datetime]Project_simulation_metadata.json'
      % 3: Other outputs
      out_name = sprintf('%s_ct%d_obj%d_data%d', simulation, jj, ii, data_index)
      disp('Results will be saved in:')
      disp(out_name)
      
      if strcmp(simulation, 'DFA')==1, 
        disp('DFA...')
        %rho = 0
        % change culture medium KSOM_AA_Jin or DMEMF12_MGSA
        DFAinterface(model_path, obj, obj_type, obj_c, save_root_path, data_path, prefix, early_stage, late_stage, out_name, ctrl, kappa, genekoflag, rxnkoflag, medium, media_perturbation, FVAflag, model_name);
      else,
      %simulation=='CFR' | simulation=='model',
        % CFR with fluxes
        if input_objective_weights==0,
          CFRinterface(model_path, pfba, obj, obj_type, obj_c, save_root_path, data_path, out_name, late_stage, early_stage, ctrl, kappa, rho, medium, genekoflag, rxnkoflag, media_perturbation, FVAflag, model_name);
        else,
          CFRinterface(model_path, pfba, obj, obj_type, obj_c, save_root_path, data_path, out_name, get_late_stage, early_stage, ctrl, kappa, rho, medium, genekoflag, rxnkoflag, media_perturbation, FVAflag, model_name);
        end
      end

    end
  end





