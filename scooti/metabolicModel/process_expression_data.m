function [ups_table, dws_table] = process_expression_data(ups, dws, model_name)
  % Choose BiGG mapping path based on model, resolving paths relative to this file
  this_dir = fileparts(mfilename('fullpath'));
  gems_dir = fullfile(this_dir, 'GEMs');
  switch model_name
    case 'Recon3D'
      map_path = fullfile(gems_dir, 'Recon3D_genes.json');
    case 'Recon2.2'
      map_path = fullfile(gems_dir, 'Recon2.2_symbol_to_hgnc.json');
    case 'Recon1'
      % Prefer user-provided Recon1 genes JSON (Recon3D-style) if present
      cand_genes = fullfile(gems_dir, 'Recon1_genes.json');
      cand_map   = fullfile(gems_dir, 'Recon1_symbol_to_hgnc.json');
      if exist(cand_genes, 'file')
        map_path = cand_genes;
      elseif exist(cand_map, 'file')
        map_path = cand_map;
      else
        warning('Recon1 mapping JSON not found; using uppercase gene symbols without conversion.');
        ups_table = ups; dws_table = dws; return;
      end
    otherwise
      warning('Model name not recognized for BiGG mapping. Skipping annotation.');
      ups_table = ups; dws_table = dws; return;
  end

  % Extract gene symbol columns robustly and uppercase
  ups_syms = local_extract_symbols(ups);
  dws_syms = local_extract_symbols(dws);

  % Convert using mapping (symbol -> BiGG ID used by model)
  up_genes = BiGG_to_annotation(ups_syms, map_path, 'to_BiGG');
  dw_genes = BiGG_to_annotation(dws_syms, map_path, 'to_BiGG');

  % Wrap back into tables
  ups_table = array2table([(1:length(up_genes))', up_genes], ...
      'VariableNames', {'Index', 'Gene'});
  dws_table = array2table([(1:length(dw_genes))', dw_genes], ...
      'VariableNames', {'Index', 'Gene'});

end

function syms = local_extract_symbols(T)
  % Return an uppercase string array of symbols from a table that may have
  % 1 or more columns; prefer a column named 'upgenes'/'dwgenes' or take first.
  if istable(T)
    vars = T.Properties.VariableNames;
    pick = 1;
    for v = 1:numel(vars)
      n = lower(string(vars{v}));
      if contains(n, 'upgenes') || contains(n, 'dwgenes') || contains(n, 'gene')
        pick = v; break;
      end
    end
    syms = upper(string(T{:,pick}));
  elseif iscell(T)
    syms = upper(string(T(:)));
  else
    syms = upper(string(T));
  end
  syms = syms(:);
end
