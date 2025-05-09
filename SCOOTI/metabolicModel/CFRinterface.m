addpath(genpath('utils'));

function CFRinterface(model_path, pfba, obj, obj_type, obj_c, root_path, data_path, out_name, upsheet, dwsheet, ctrl, kappa, rho, medium, genekoflag, rxnkoflag, media_perturbation, FSflag, model_name, CFR_model, extra_weight, algorithm)

  sample_name = sprintf('%s%s', upsheet, dwsheet);

%% Load and config metabolic network model
  model = load_and_configure_model(model_path, model_name, medium);

%% manage file names
  % File name with datetime as prefix
  file_prefix = string(datetime('now','TimeZone','local','Format','MMMdyHHmm'));
  % Output file
  filename = sprintf('%s/[%s]%s', root_path, file_prefix, out_name);
  disp(sprintf('%s %s', 'The CFR result has been saved in', filename));
  % Output file
  excelname = filename;
  
  
%% Set objectives for multi- or single-objective problems
  model = objective_setting_function(model, obj, obj_c, obj_type, model_name, algorithm);
  save_metadata_json(excelname, model, obj, obj_type, obj_c, ...
    root_path, data_path, out_name, ctrl, kappa, rho, medium, ...
    genekoflag, rxnkoflag, media_perturbation, upsheet, dwsheet, ...
    model_path, CFR_model, extra_weight, algorithm);
 
  
  % Separate values and labels
  if ctrl==0,
    ups = {};
    dws = {};
  else,
    [ups, dws] = load_expression_tables(data_path, upsheet, dwsheet);
  end
 
%% manage model parameters
  size_ups = size(ups); size_dws = size(dws);
  size_ups = size_ups(1); size_dws = size_dws(1); 
  if size_ups>0 & size_dws>0,
    [ups, dws] = process_expression_data(ups, dws, model_name);
  end

  %% main function for CFR
 run_simulation(model, excelname, algorithm, ups, dws, ...
               ctrl, FSflag, kappa, rho, ...
               recon_model, extra_weight, ...
               genekoflag, rxnkoflag, file_prefix, ...
               out_name, root_path, metadata);
  % end if run simulation for CFR 
end
% end if CFRinterface
