function model = load_config_model(model_path, model_name, medium)
  % Load the metabolic network model
  raw = load(model_path);
  fn = fieldnames(raw);
  model = raw.(fn{1});

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
    for i = 1:length(EX_rxns)
      rxnID = findRxnIDs(model, EX_rxns{i});
      if rxnID > 0
        model = changeRxnBounds(model, model.rxns{rxnID}, EX_mets{i}, 'l');
      end
    end
  end
end

