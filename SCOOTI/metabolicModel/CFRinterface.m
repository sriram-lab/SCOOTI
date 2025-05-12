function CFRinterface(config)
% CFRinterface - Main entry for constraint-based flux regulation using config struct input

  addpath(genpath('utils'));
  sample_name = sprintf('%s%s', config.uplabel, config.dwlabel);

  %% Load and configure metabolic network model
  model = load_config_model(config.GEM_path, config.model_name, config.medium);

  %% Manage file names
  file_prefix = string(datetime('now','TimeZone','local','Format','MMMdyHHmm'));
  filename = sprintf('%s/[%s]%s', config.save_root_path, file_prefix, config.out_name);
  disp(sprintf('%s %s', 'The CFR result has been saved in', filename));
  excelname = filename;


  %% Set objectives
  model = objective_setting_function(model, config.obj, config.obj_c, config.obj_type, config.model_name, config.algorithm);
  config.metadata = save_metadata_json(excelname, model, config);
  %config.obj, config.obj_type, config.obj_c, ...
  %  config.save_root_path, config.data_path, config.out_name, config.ctrl, config.kappa, config.rho, config.medium, ...
  %  config.genekoflag, config.rxnkoflag, config.medium_perturbation, config.uplabel, config.dwlabel, ...
  %  config.GEM_path, config.CFR_model_path, config.extraWeight, config.algorithm);

  %% Load expression data
  if config.ctrl == 0
    ups = {};
    dws = {};
  else
    [ups, dws] = load_expression_tables(config.data_path, config.uplabel, config.dwlabel);
  end

  %% load CFR model
  if ~isempty(config.CFR_model_path)
    config.CFR_models = load_CFR_models(CFR_model_path)
  end

  %% Manage model parameters
  if ~isempty(ups) && ~isempty(dws)
    [ups, dws] = process_expression_data(ups, dws, config.model_name);
    disp(ups)
  end

  %% Main simulation run
  run_simulation(model, excelname, config.algorithm, ups, dws, ...
                 config.ctrl, config.FSflag, config.kappa, config.rho, ...
                 config.CFR_model, config.extraWeight, ...
                 config.genekoflag, config.rxnkoflag, file_prefix, ...
                 config.out_name, config.save_root_path, config.metadata);
end




