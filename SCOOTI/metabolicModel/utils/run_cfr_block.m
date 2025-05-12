function run_cfr_block(model, uplist, dwlist, excelname, ...
                       sample_name, out_name, root_path, ...
                       kappa, rho, genekoflag, rxnkoflag, ...
                       recon_model, extra_weight)

  % Run CFR with constraints
  [uncon_flux, fluxstate, grate_naive, geneko_flux, rxnko_growthrate, ...
   solverobj_naive, model_out] = constrain_flux_regulation( ...
      model, uplist, dwlist, kappa, rho, 1E-3, 0, ...
      genekoflag, rxnkoflag, [], [], recon_model, extra_weight);

  % Save context-specific model
  model_filename = sprintf('%s_%s.mat', excelname, 'model_CFR');
  save(model_filename, "model_out");

  % Choose flux to write based on pfba flag (default to fluxstate)
  gh = fluxstate;
  gh(end+1, 1) = grate_naive;

  % Prepare output table
  rxns = model.rxns;
  rxns{end+1} = 'Obj';
  T = table(rxns, gh, 'VariableNames', {'rxns', sample_name});

  % Save flux results
  excelsheet = 'fluxes';
  tmp_filename = sprintf('%s_%s.csv', excelname, excelsheet);
  writetable(T, tmp_filename);
  fprintf('[CFR] Flux results saved to %s\n', tmp_filename);

  % Gene KO analysis
  if genekoflag
    unigenes = unique(model.genes);
    unigenes(cellfun(@isempty, unigenes)) = {'unknown'};
    genes = ['WT', unigenes(:)'];
    gene_deletions = [model.rxns, num2cell(geneko_flux)];
    mat_filename = sprintf('%s_%s.mat', excelname, 'CFR-geneDel');
    save(mat_filename, "rxns", "genes", "geneko_flux");
    fprintf('[CFR] Gene deletion results saved to %s\n', mat_filename);
  end

  % Reaction KO logic can be added here if needed
end

