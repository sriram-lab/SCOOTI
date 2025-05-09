function run_simulation(model, excelname, algorithm, ups, dws, ...
                        ctrl, FSflag, kappa, rho, ...
                        recon_model, extra_weight, ...
                        genekoflag, rxnkoflag, file_prefix, ...
                        out_name, root_path, metadata)

  % Initialize result table
  T = table();
  rxns = model.rxns;
  rxns{end+1} = 'Obj';  % Append objective row
  sample_name = 'var';  % General placeholder

  % Determine up/down gene list based on model
  genes = model.genes;
  uk = unique(genes);

  if ctrl == 0 || strcmp(algorithm, 'MOOMIN')
    uplist = {};
    dwlist = {};
  else
    uplist = filter_gene_list(ups, uk);
    dwlist = filter_gene_list(dws, uk);
  end

  if FSflag == 1
    % Run flux sampling (e.g., CHRR)
    run_sampling_block(model, excelname, root_path, file_prefix, ...
                       out_name, metadata);
  elseif strcmp(algorithm, 'INIT')
    run_init_block(model, uplist, dwlist, excelname, file_prefix, ...
                   out_name, root_path, metadata);
  elseif strcmp(algorithm, 'MOOMIN')
    run_moomin_block(model, ups, excelname, sample_name, ...
                     out_name, root_path, genekoflag, metadata);
  else
    run_cfr_block(model, uplist, dwlist, excelname, ...
                  sample_name, out_name, root_path, ...
                  kappa, rho, genekoflag, rxnkoflag, ...
                  recon_model, extra_weight, metadata);
  end
end

