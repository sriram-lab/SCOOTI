function config = set_default_parameters(config)
  % Fill missing fields in config with sensible defaults



  % Backward compatibility: map generic kappa/rho keys to CFR-specific keys
  % if provided in JSON configs.
  if (~isfield(config, 'CFR_kappa') || isempty(config.CFR_kappa)) && isfield(config, 'kappa') && ~isempty(config.kappa)
    config.CFR_kappa = config.kappa;
  end
  if (~isfield(config, 'CFR_rho') || isempty(config.CFR_rho)) && isfield(config, 'rho') && ~isempty(config.rho)
    config.CFR_rho = config.rho;
  end

  if ~isfield(config, 'jj') || isempty(config.jj), config.jj = 1; end
  if ~isfield(config, 'input_obj_tb') || isempty(config.input_obj_tb), config.input_obj_tb = ''; end
  if ~isfield(config, 'DFA_kappa') || isempty(config.DFA_kappa), config.DFA_kappa = 1; end
  if ~isfield(config, 'CFR_kappa') || isempty(config.CFR_kappa), config.CFR_kappa = 1; end
  if ~isfield(config, 'CFR_rho') || isempty(config.CFR_rho), config.CFR_rho = 1; end
  if ~isfield(config, 'extraWeight') || isempty(config.extraWeight), config.extraWeight = 1; end
  if ~isfield(config, 'pairwise_CFR_model') || isempty(config.pairwise_CFR_model), config.pairwise_CFR_model = 0; end
  if ~isfield(config, 'init_objective') || isempty(config.init_objective), config.init_objective = 1; end
  if ~isfield(config, 'algorithm') || isempty(config.algorithm), config.algorithm = 'iMAT'; end
  if ~isfield(config, 'genekoflag') || isempty(config.genekoflag), config.genekoflag = 0; end
  if ~isfield(config, 'rxnkoflag') || isempty(config.rxnkoflag), config.rxnkoflag = 0; end
  if ~isfield(config, 'FSflag') || isempty(config.FSflag), config.FSflag = 0; end
  if ~isfield(config, 'medium_perturbation') || isempty(config.medium_perturbation), config.medium_perturbation = 0; end
  if ~isfield(config, 'pfba') || isempty(config.pfba), config.pfba = 1; end
  if ~isfield(config, 'late_stage') || isempty(config.late_stage), config.late_stage = 'upgenes'; end
  if ~isfield(config, 'early_stage') || isempty(config.early_stage), config.early_stage = 'dwgenes'; end
end



%function config = set_default_parameters(config)
%  defaults = struct( ...
%    'jj', 1, 'input_obj_tb', '', 'DFA_kappa', 1, 'CFR_kappa', 1, 'CFR_rho', 1, ...
%    'extraWeight', 1, 'pairwise_CFR_model', 0, 'init_objective', 1, 'algorithm', 'iMAT', ...
%    'genekoflag', 0, 'rxnkoflag', 0, 'FSflag', 0, 'medium_perturbation', 0, 'pfba', 1, ...
%    'late_stage', 'upgenes', 'early_stage', 'dwgenes', 'input_objective_weights', 0);
%
%  fields = fieldnames(defaults);
%  for i = 1:numel(fields)
%    if ~isfield(config, fields{i}) || isempty(config.(fields{i}))
%      config.(fields{i}) = defaults.(fields{i});
%    end
%  end
