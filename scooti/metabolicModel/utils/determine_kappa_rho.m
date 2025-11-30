function [kappa, rho] = determine_kappa_rho(config)
% determine_kappa_rho - decide kappa/rho for DFA/CFR based on config
% Supports:
%  - Fixed numeric values (CFR_kappa, CFR_rho)
%  - Logspace scanning when *_kappa == -1 using paraLen/random_para/jj
%  - Explicit lists via strings or arrays:
%      kappaArr/rhoArr as comma-separated strings or numeric vectors
%      CFR_kappa/CFR_rho as comma-separated strings

  %% Helpers to coerce possible list fields to numeric row vectors
  function arr = parse_list_field(val)
    if isstring(val) || ischar(val)
      if contains(string(val), ",")
        parts = split(string(val), ",");
        arr = str2double(parts(:))';
      else
        % single numeric encoded as string
        arr = str2double(string(val));
      end
    elseif isnumeric(val)
      if isscalar(val)
        arr = val;
      else
        arr = val(:)';
      end
    elseif iscell(val)
      arr = cellfun(@(x) double(string(x)), val);
      arr = arr(:)';
    else
      arr = [];
    end
  end

  %% DFA branch
  if strcmp(config.simulation, 'DFA')
    % Optional list via dkappaArr
    k_list = [];
    if isfield(config, 'dkappaArr')
      k_list = parse_list_field(config.dkappaArr);
    end
    if isempty(k_list) && (isfield(config, 'DFA_kappa') && (isstring(config.DFA_kappa) || ischar(config.DFA_kappa)))
      k_list = parse_list_field(config.DFA_kappa);
    end

    if ~isempty(k_list)
      % Use provided list
      N = numel(k_list);
      sel = min(max(1, config.jj), N);
      if isfield(config,'random_para') && config.random_para
        sel = randi(N);
      end
      kappa = k_list(sel);
      rho = 0;
      return;
    end

    if isfield(config,'DFA_kappa') && config.DFA_kappa == -1
      space = logspace(1, -3, config.paraLen);
      if isfield(config,'random_para') && config.random_para
        idx = randsample(numel(space), min(config.paraLen, numel(space)));
        space = space(idx);
      end
      sel = min(max(1, config.jj), numel(space));
      kappa = space(sel);
      rho = 0;
    else
      kappa = config.DFA_kappa;
      rho = 0;
    end
    return;
  end

  %% CFR branch (or model)
  % Accept explicit lists from kappaArr/rhoArr or comma-separated CFR_kappa/CFR_rho
  k_list = [];
  r_list = [];
  if isfield(config, 'kappaArr')
    k_list = parse_list_field(config.kappaArr);
  end
  if isfield(config, 'rhoArr')
    r_list = parse_list_field(config.rhoArr);
  end
  if isempty(k_list) && isfield(config, 'CFR_kappa') && (isstring(config.CFR_kappa) || ischar(config.CFR_kappa))
    k_list = parse_list_field(config.CFR_kappa);
  end
  if isempty(r_list) && isfield(config, 'CFR_rho') && (isstring(config.CFR_rho) || ischar(config.CFR_rho))
    r_list = parse_list_field(config.CFR_rho);
  end

  if ~isempty(k_list) || ~isempty(r_list)
    % Build parameter grid according to provided lists; support fixed-other
    if isempty(k_list)
      k_list = config.CFR_kappa;
    end
    if isempty(r_list)
      r_list = config.CFR_rho;
    end
    [Y, Z] = meshgrid(k_list, r_list);
    param_grid = [Y(:), Z(:)];
    % Optionally subsample to paraLen items
    total = size(param_grid,1);
    if isfield(config,'paraLen') && config.paraLen > 0 && config.paraLen < total
      if isfield(config,'random_para') && config.random_para
        idx = randsample(total, config.paraLen);
      else
        idx = 1:config.paraLen;
      end
      param_grid = param_grid(idx, :);
    end
    sel = min(max(1, config.jj), size(param_grid,1));
    kappa = param_grid(sel,1);
    rho   = param_grid(sel,2);
    return;
  end

  % Legacy: logspace grid when CFR_kappa == -1
  if isfield(config,'CFR_kappa') && config.CFR_kappa == -1
    space = logspace(1, -3, config.paraLen);
    [Y, Z] = meshgrid(space, space);
    param_grid = [Y(:), Z(:)];
    if isfield(config,'random_para') && config.random_para
      idx = randsample(size(param_grid, 1), min(config.paraLen, size(param_grid,1)));
      param_grid = param_grid(idx, :);
    end
    sel = min(max(1, config.jj), size(param_grid,1));
    kappa = param_grid(sel, 1);
    rho   = param_grid(sel, 2);
  else
    kappa = config.CFR_kappa;
    rho   = config.CFR_rho;
  end
end
