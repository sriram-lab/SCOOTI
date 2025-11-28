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
    % Extract raw gene lists (uppercase) for diagnostics
    raw_up = local_extract_genes(ups);
    raw_dw = local_extract_genes(dws);
    % Filter against model genes (already uppercased upstream)
    uplist = filter_gene_list(ups, uk);
    dwlist = filter_gene_list(dws, uk);

    % Diagnostics: matching statistics
    try
      nu = numel(raw_up); nd = numel(raw_dw);
      mu = numel(uplist); md = numel(dwlist);
      fprintf('[CFR] Gene matching (UP): %d input -> %d matched in model\n', nu, mu);
      fprintf('[CFR] Gene matching (DW): %d input -> %d matched in model\n', nd, md);
      if mu == 0
        % Show a few examples that failed to match
        up_miss = setdiff(unique(raw_up), unique(uplist));
        k = min(10, numel(up_miss));
        if k > 0
          fprintf('[CFR][warn] Example UP genes not in model (first %d):\n', k);
          disp(up_miss(1:k));
        end
      end
      if md == 0
        dw_miss = setdiff(unique(raw_dw), unique(dwlist));
        k = min(10, numel(dw_miss));
        if k > 0
          fprintf('[CFR][warn] Example DW genes not in model (first %d):\n', k);
          disp(dw_miss(1:k));
        end
      end
    catch
      % ignore diagnostics errors
end


function genes = local_extract_genes(gene_table)
  % Extract an uppercase gene list from a table/cell/array without intersecting with model
  genes = strings(0,1);
  try
    if istable(gene_table)
      w = width(gene_table);
      if w >= 2
        genes = upper(string(gene_table{:,2}));
      elseif w >= 1
        genes = upper(string(gene_table{:,1}));
      else
        genes = strings(0,1);
      end
    elseif iscell(gene_table)
      genes = upper(string(gene_table(:)));
    else
      genes = upper(string(gene_table));
    end
  catch
    genes = strings(0,1);
  end
  genes = genes(:);
end

    % Convert uplist safely
    if isstring(uplist)
        uplist = cellstr(uplist);
    end
    
    % Also ensure dwlist
    if isstring(dwlist)
        dwlist = cellstr(dwlist);
    end


    % Preview a few matched UP genes
    try
      k = min(10, numel(uplist));
      if k>0
        fprintf('[CFR] Preview matched UP genes (first %d):\n', k);
        disp(uplist(1:k));
      else
        fprintf('[CFR][warn] No matched UP genes found. Constraints may be ineffective.\n');
      end
    catch
    end
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
                     out_name, root_path, genekoflag);
  else
    run_cfr_block(model, uplist, dwlist, excelname, ...
                  sample_name, out_name, root_path, ...
                  kappa, rho, genekoflag, rxnkoflag, ...
                  recon_model, extra_weight);
  end
end
