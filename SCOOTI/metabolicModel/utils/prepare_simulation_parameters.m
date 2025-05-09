function params = prepare_simulation_parameters(config)
  idx = config.data_idx;

  params.model = config.model;
  params.data_path = config.data_series{idx};
  params.medium = config.medium_series{idx};
  params.prefix = config.prefix_series{idx};
  params.ctrl = config.ctrl;
  params.kappa = config.kappa;
  params.rho = config.rho;
  params.FSflag = config.FSflag;
  params.genekoflag = config.genekoflag;
  params.rxnkoflag = config.rxnkoflag;
  params.medium_perturbation = config.medium_perturbation;
  params.save_root_path = config.save_root_path;
  params.model_name = config.model_name;
  params.pfba = config.pfba;
  params.extraWeight = config.extraWeight;
  params.algorithm = config.algorithm;

  % Objective setup
  if config.input_objective_weights == 0
    ii = config.init_objective;
  else
    ii = 2; % default to column 2 if multi-objective input
  end

  [params.obj, params.obj_c] = select_objectives(config, ii);

  % Output name convention
  params.out_name = sprintf('%s_ct%d_data%d', config.simulation, config.jj, idx);
end




%function params = prepare_simulation_parameters(config)
%  idx = config.data_idx;
%
%  params.model = config.model;
%  params.data_path = config.data_series{idx};
%  params.medium = config.medium_series{idx};
%  params.prefix = config.prefix_series{idx};
%  params.ctrl = config.ctrl;
%  params.kappa = config.kappa;
%  params.rho = config.rho;
%  params.FSflag = config.FSflag;
%  params.genekoflag = config.genekoflag;
%  params.rxnkoflag = config.rxnkoflag;
%  params.medium_perturbation = config.medium_perturbation;
%  params.save_root_path = config.save_root_path;
%  params.model_name = config.model_name;
%  params.pfba = config.pfba;
%  params.extraWeight = config.extraWeight;
%  params.algorithm = config.algorithm;
%
%  % Output name convention
%  params.out_name = sprintf('%s_ct%d_data%d', config.simulation, config.jj, idx);
%end


