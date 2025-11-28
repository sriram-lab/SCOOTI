function CFRinterface(config)
% CFRinterface - Main entry for constraint-based flux regulation using config struct input

  addpath(genpath('utils'));
  % Resolve repository root from this file's location so we can
  % build absolute output paths regardless of current directory
  this_dir = fileparts(mfilename('fullpath'));
  repo_root = fileparts(fileparts(this_dir)); % .../scooti/metabolicModel -> repo root
  sample_name = sprintf('%s%s', config.uplabel, config.dwlabel);

  %% Load and configure metabolic network model
  model = load_config_model(config.GEM_path, config.model_name, config.medium);
  % (Reverted additions: avoid altering model containers or printing model genes here)

  %% Manage file names
  file_prefix = string(datetime('now','TimeZone','local','Format','MMMdyHHmm'));
  save_dir = fullfile(repo_root, config.save_root_path);
  if ~exist(save_dir, 'dir')
    mkdir(save_dir);
  end
  filename = sprintf('%s/[%s]%s', save_dir, file_prefix, config.out_name);
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
    % Print a compact preview to avoid flooding the console
    try
      fprintf('[CFRinterface] Preview of upregulated genes (first 10 rows):\n');
      disp(head_safe(ups, 10));
      fprintf('[CFRinterface] Preview of downregulated genes (first 10 rows):\n');
      disp(head_safe(dws, 10));
    catch
      % Fallback: print sizes only
      try
        fprintf('[CFRinterface] ups size: %dx%d\n', size(ups,1), size(ups,2));
        fprintf('[CFRinterface] dws size: %dx%d\n', size(dws,1), size(dws,2));
      catch
        % ignore if size fails
      end
    end
  end

  %% Main simulation run
  run_simulation(model, excelname, config.algorithm, ups, dws, ...
                 config.ctrl, config.FSflag, config.kappa, config.rho, ...
                 config.CFR_model, config.extraWeight, ...
                 config.genekoflag, config.rxnkoflag, file_prefix, ...
                 config.out_name, config.save_root_path, config.metadata);
end




function out = head_safe(x, n)
% head_safe - return first n rows of a table/cell/array without errors
  if nargin < 2, n = 10; end
  out = x;
  try
    if istable(x)
      n = min(n, height(x));
      out = x(1:n, :);
    elseif iscell(x)
      n = min(n, size(x,1));
      out = x(1:n, :);
    elseif isnumeric(x) || islogical(x)
      n = min(n, size(x,1));
      out = x(1:n, :);
    else
      % Unknown type; return as-is
      out = x;
    end
  catch
    out = x;
  end
end
