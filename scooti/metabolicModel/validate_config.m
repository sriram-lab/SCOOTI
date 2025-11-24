function validate_config(config)
% validate_config - Ensures required fields in config are present and valid
% Throws error if any required field is missing or malformed.

  required_fields = {
    'GEM_path', 'COBRA_path', 'model_name', ...
    'save_root_path', 'simulation', 'prefix_name'
  };

  for i = 1:numel(required_fields)
    field = required_fields{i};
    if ~isfield(config, field) || isempty(config.(field))
      error('Missing or empty required config field: %s', field);
    end
  end

  % Optional: check file existence
  if ~exist(config.GEM_path, 'file')
    error('GEM_path does not exist: %s', config.GEM_path);
  end
  if ~exist(config.COBRA_path, 'dir')
    error('COBRA_path does not exist: %s', config.COBRA_path);
  end

  % Optional: check types
  if ~ischar(config.model_name) && ~isstring(config.model_name)
    error('model_name must be a string');
  end

  % Optional: issue warning for optional fields
  if ~isfield(config, 'obj_candidate_list_file') || isempty(config.obj_candidate_list_file)
    warning('obj_candidate_list_file is empty; default objectives may be used');
  end
end

