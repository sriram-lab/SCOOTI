function CFR_models = load_CFR_models(CFR_model_path)
  CFR_models = {};
  if ~isempty(CFR_model_path)
    files = dir(CFR_model_path);
    for i = 1:length(files)
      if contains(files(i).name, '.mat')
        CFR_models{end+1, 1} = fullfile(CFR_model_path, files(i).name);
      end
    end
  end
end
