function run_moomin_block(model, ups, excelname, sample_name, ...
                          out_name, root_path, genekoflag, metadata)

  % Fix model rules before MOOMIN
  model.b = model.b(:, 1);
  model = fixModelRules(model);

  rxns = model.rxns;
  rxns{end+1} = 'Obj';

  if ~genekoflag
    % Just run MOOMIN once for WT
    [fluxstate, objFlux] = run_moomin(model, ups);
    gh = fluxstate;
    gh(end+1) = objFlux;

    T = table(rxns, gh, 'VariableNames', {'rxns', sample_name});
    tmp_filename = sprintf('%s_fluxes.csv', excelname);
    writetable(T, tmp_filename);
    fprintf('[MOOMIN] Flux result saved to %s\n', tmp_filename);

  else
    % Run gene KO analysis with MOOMIN
    fprintf('[MOOMIN] Running gene knockouts...\n');
    unqgenes = unique(model.genes);
    save_genes = [{'WT'}, unqgenes(:)'];
    geneko_flux = zeros(length(model.rxns), length(save_genes));

    % WT fluxes
    [fluxstate, objFlux] = run_moomin(model, ups);
    geneko_flux(:, 1) = fluxstate;

    % KO fluxes
    for k = 1:length(unqgenes)
      gene = unqgenes{k};
      try
        ko_model = deleteModelGenes(model, gene);
        ko_model.b = ko_model.b(:, 1);
        ko_model = fixModelRules(ko_model);
        fluxstate = run_moomin(ko_model, ups);
        geneko_flux(:, k+1) = fluxstate;
      catch
        warning('Could not KO gene: %s', gene);
        geneko_flux(:, k+1) = NaN;
      end
    end

    % Clean unknown gene names
    save_genes(cellfun(@isempty, save_genes)) = {'unknown'};

    % Save as .mat
    mat_filename = sprintf('%s_CFR-geneDel.mat', excelname);
    save(mat_filename, "rxns", "save_genes", "geneko_flux");
    fprintf('[MOOMIN] Gene deletion results saved to %s\n', mat_filename);
  end
end

