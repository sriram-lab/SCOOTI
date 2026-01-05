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
  M = local_load_mapping(map_path);
  up_genes = strings(numel(ups_syms),1);
  dw_genes = strings(numel(dws_syms),1);
  for i=1:numel(ups_syms)
    key = upper(strtrim(ups_syms(i)));
    if isKey(M, key); up_genes(i) = string(M(key)); else; up_genes(i) = string(ups_syms(i)); end
  end
  for i=1:numel(dws_syms)
    key = upper(strtrim(dws_syms(i)));
    if isKey(M, key); dw_genes(i) = string(M(key)); else; dw_genes(i) = string(dws_syms(i)); end
  end
  up_genes = cellstr(up_genes);
  dw_genes = cellstr(dw_genes);

  % Wrap back into tables
  ups_table = table((1:numel(up_genes))', up_genes, 'VariableNames', {'Index','Gene'});
  dws_table = table((1:numel(dw_genes))', dw_genes, 'VariableNames', {'Index','Gene'});

end

function M = local_load_mapping(p)
  % Return containers.Map mapping SYMBOL (upper) -> ID (as-is)
  fid = fopen(p, 'r');
  if fid < 0, error('Cannot open mapping JSON: %s', p); end
  raw = fread(fid, inf, '*char')'; fclose(fid);
  % Try structured JSON first
  try
    data = jsondecode(raw);
  catch
    data = [];
  end
  names = [];
  ids = [];
  if ~isempty(data)
    if isstruct(data) && isfield(data, 'results')
      names = {data.results.name};
      ids   = {data.results.bigg_id};
    elseif isstruct(data) && isfield(data, 'symbol') && isfield(data, 'hgnc')
      names = data.symbol;
      ids   = data.hgnc;
    end
  end
  if isempty(names) || isempty(ids)
    % Fallback: parse flat object mapping {"BiGG_id":"SYMBOL", ...}
    toks = regexp(raw, '"([^"\\]+)"\s*:\s*"([^"\\]*)"', 'tokens');
    if isempty(toks)
      error('Unsupported mapping schema in %s', p);
    end
    % Invert: symbol -> BiGG_id
    ids = cellfun(@(t) t{1}, toks, 'UniformOutput', false);      % keys (BiGG ids)
    names = cellfun(@(t) t{2}, toks, 'UniformOutput', false);    % values (symbols)
  end
  names = upper(string(names(:)));
  ids   = string(ids(:));
  M = containers.Map('KeyType','char','ValueType','char');
  for i=1:numel(names)
    k = char(names(i)); v = char(ids(i));
    if ~isempty(k) && ~isempty(v) && ~isKey(M, k)
      M(k) = v;
    end
  end
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
