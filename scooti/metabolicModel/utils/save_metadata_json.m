function metadata = save_metadata_json(filename, model, config)
%function metadata = save_metadata_json(filename, model, obj, obj_type, obj_c, ...
%    root_path, data_path, out_name, ctrl, kappa, rho, medium, ...
%    genekoflag, rxnkoflag, media_perturbation, upsheet, dwsheet, ...
%    model_path, CFR_model, extra_weight, algorithm)


  %disp(config)

  % Extract nonzero objectives
  obj_indices = find(model.c ~= 0);
  if ~isempty(obj_indices)
    obj_weights = model.c(obj_indices);
    obj_rxns = {model.rxns(obj_indices)};
  else
    obj_weights = 0;
    obj_rxns = '';
  end

  % Decide sample name
  if config.ctrl == 1
    input_path = config.data_path;
  else
    input_path = sprintf('%s%s', config.uplabel, config.dwlabel);
  end

  % Build metadata struct
  metadata.obj = config.obj
  metadata.obj_type = config.obj_type
  metadata.obj_c = config.obj_c
  metadata.output_path = config.save_root_path
  metadata.input_path = input_path
  metadata.file_name = config.out_name
  metadata.with_constraint = config.ctrl
  metadata.CFR_kappa = config.kappa
  metadata.CFR_rho = config.rho
  metadata.medium = config.medium
  metadata.genekoflag = config.genekoflag
  metadata.rxnkoflag = config.rxnkoflag
  metadata.media_perturbation = config.medium_perturbation
  metadata.objWeights = obj_weights
  metadata.objRxns = obj_rxns
  metadata.model_path = config.GEM_path
  metadata.upStage = config.uplabel
  metadata.dwStage = config.dwlabel
  metadata.CFRModel = config.CFR_model_path
  metadata.extraWeight = config.extraWeight
  metadata.algorithm = config.algorithm

  % Encode and write to JSON file
  json_str = jsonencode(metadata);
  json_file = sprintf('%s_metadata.json', filename);
  fid = fopen(json_file, 'w');
  fprintf(fid, '%s', json_str);
  fclose(fid);
end

