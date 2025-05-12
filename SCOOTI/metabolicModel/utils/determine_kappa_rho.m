function [kappa, rho] = determine_kappa_rho(config)
  if strcmp(config.simulation, 'DFA')
    if config.DFA_kappa == -1
      space = logspace(1, -3, config.paraLen);
      if config.random_para
        idx = randsample(numel(space), config.paraLen);
        space = space(idx);
      end
      kappa = space(config.jj);
      rho = 0;
    else
      kappa = config.DFA_kappa;
      rho = 0;
    end
  else % CFR or model
    if config.CFR_kappa == -1
      space = logspace(1, -3, config.paraLen);
      [Y, Z] = meshgrid(space, space);
      param_grid = [Y(:), Z(:)];
      if config.random_para
        idx = randsample(size(param_grid, 1), config.paraLen);
        param_grid = param_grid(idx, :);
      end
      kappa = param_grid(config.jj, 1);
      rho = param_grid(config.jj, 2);
    else
      kappa = config.CFR_kappa;
      rho = config.CFR_rho;
    end
  end
end
