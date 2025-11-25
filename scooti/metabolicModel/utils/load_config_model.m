function model = load_config_model(model_path, model_name, medium)
  % Load the metabolic network model
  raw = load(model_path);
  fn = fieldnames(raw);
  model = raw.(fn{1});

  % Normalize model fields to COBRA-expected types (cell arrays of char)
  normFields = {'rxns','rxnNames','mets','metNames','genes','grRules','subSystems'};
  for k = 1:numel(normFields)
    f = normFields{k};
    if isfield(model, f)
      v = model.(f);
      try
        if isstring(v)
          model.(f) = cellstr(v);
        elseif iscell(v)
          % Convert any string scalars to char within the cell
          model.(f) = cellfun(@(x) char(x), v, 'UniformOutput', false);
        elseif ischar(v)
          model.(f) = {v};
        end
      catch
        % Best-effort fallback via string() then cellstr
        try
          model.(f) = cellstr(string(v));
        catch
          % leave as-is
        end
      end
    end
  end


  % Special handling for Recon1
  if strcmp(model_name, 'Recon1')
    bio_obj_idx = find(contains(model.rxns, 'biomass_objective'));
    if ~isempty(bio_obj_idx)
      model.c(bio_obj_idx) = 1;
    end
  end

  % Modify medium if specified
  if ~isempty(medium)
    media = readtable('FINAL_MEDIUM_MAP_RECON1.xlsx', 'Sheet', medium);
    EX_rxns = media{:,6};
    EX_mets = media{:,3};
    % Coerce bounds to numeric (Excel may load as cell array of strings)
    if iscell(EX_mets)
      tmp = nan(size(EX_mets));
      for k = 1:numel(EX_mets)
        v = EX_mets{k};
        if ischar(v) || isstring(v)
          vn = str2double(string(v));
        elseif isnumeric(v)
          vn = v;
        else
          vn = NaN;
        end
        if isnan(vn); vn = 0; end
        tmp(k) = vn;
      end
      EX_mets = tmp;
    elseif ~isnumeric(EX_mets)
      try
        EX_mets = double(EX_mets);
      catch
        % As a fallback, set to zeros with warning
        warning('EX_mets column not numeric; defaulting to zeros.');
        EX_mets = zeros(size(EX_rxns));
      end
    end
    for i = 1:length(EX_rxns)
      rxnID = findRxnIDs(model, EX_rxns{i});
      if rxnID > 0
        model = changeRxnBounds(model, model.rxns{rxnID}, EX_mets(i), 'l');
      end
    end
  end
end
