function test_run_CFR_pipeline()
  config = struct();
  config.model_path = 'models/Recon1.mat';
  config.obj = {'biomass_objective'};
  config.obj_type = {'min'};
  config.obj_c = [1];
  config.root_path = 'results';
  config.data_path = 'expression_data.xlsx';
  config.out_name = 'demo_run';
  config.upsheet = 'UP';
  config.dwsheet = 'DW';
  config.ctrl = true;
  config.kappa = 0.1;
  config.rho = 1;
  config.medium = "KSOM_AA";
  config.genekoflag = false;
  config.rxnkoflag = false;
  config.media_perturbation = false;
  config.FSflag = false;
  config.model_name = 'Recon1';
  config.CFR_model = '';
  config.extra_weight = 0;
  config.algorithm = 'CFR';

  run_CFR_pipeline(config);
end

