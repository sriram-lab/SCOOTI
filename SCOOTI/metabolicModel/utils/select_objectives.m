function [obj, obj_c, obj_type] = select_objectives(config, ii)
  % Decide between single-objective or multi-objective mode
  if config.input_objective_weights == 0
    % Single-objective mode
    single_obj = [{''}; config.obj_candidates{:,2}];
    obj = single_obj;
    obj_c = zeros(1, length(obj));
    obj_c(ii) = 1;
  elseif strcmp(config.prefix_name, 'model')
    % Multi-objective (unconstrained)
    obj = config.input_obj_tb{2:end, 1};
    obj_c = config.input_obj_tb{2:end, ii};
  else
    % Multi-objective (with constraint)
    obj = config.input_obj_tb{2:end, 1};
    obj_c = config.input_obj_tb{2:end, ii};
  end

  if strcmp(obj, ''),
      obj_type = '';
  else,
      obj_type = 'Demand';
  end
end


