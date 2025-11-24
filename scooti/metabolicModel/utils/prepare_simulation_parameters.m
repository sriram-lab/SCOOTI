function params = prepare_simulation_parameters(config)
  idx = config.data_idx;

  params.GEM_path = config.GEM_path;
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
  params.uplabel = config.uplabel;
  params.dwlabel = config.dwlabel;
  params.CFR_model_path = config.CFR_model_path;
  params.CFR_model = config.CFR_model;


  [params.obj, params.obj_c, params.obj_type] = select_objectives(config, config.ii);

  % Output name convention
  params.out_name = sprintf('%s_ct%d_obj%d_data%d', config.simulation, config.jj, config.ii, idx);
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


