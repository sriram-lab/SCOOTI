function run_init_block(model, uplist, dwlist, excelname, file_prefix, ...
                        out_name, root_path, metadata)

  % Prepare INIT model
  initial_model = model;
  initial_model.b = model.b(:, 1);
  initial_model = fixModelRules(initial_model);
  initial_model.c = zeros(size(initial_model.c));
  model = initial_model;

  % Run INIT to generate constrained model
  fprintf('[INIT] Running INIT algorithm...\n');
  constrained_model = run_init2(model, uplist, dwlist, 36, 1);
  [flux_avg, flux_all, rxn_inactive] = run_CHRR2(model, constrained_model, 100);

  % Output table setup
  rxns = model.rxns;
  rxns{end+1} = 'Obj';

  for col_num = 1:size(flux_all, 2)
    gh = flux_all(:, col_num);
    gh(end+1) = 0;  % Add placeholder for objective value

    T = table(rxns, gh, 'VariableNames', {'rxns', 'flux'});
    
    % Create subdirectory and filename
    subfolder = sprintf('%s/sample_%d/', root_path, col_num);
    mkdir(subfolder);
    excelname_i = sprintf('%s/sample_%d/[%s]%s', root_path, col_num, file_prefix, out_name);

    % Save metadata
    json_file = sprintf('%s_metadata.json', excelname_i);
    fid = fopen(json_file, 'w');
    fprintf(fid, '%s', jsonencode(metadata));
    fclose(fid);

    % Save table
    tmp_filename = sprintf('%s_fluxes.csv', excelname_i);
    writetable(T, tmp_filename);
    fprintf('[INIT] Flux sample %d saved to %s\n', col_num, tmp_filename);
  end
end

