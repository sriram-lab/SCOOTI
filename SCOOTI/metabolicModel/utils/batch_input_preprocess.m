function [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium)

  %% Significant genes of single cell human cell datasets
  data_dir_files = dir(data_dir);
  data_series= {};
  prefix_series = {};
  for kk=3:length(data_dir_files),
    data_name = strsplit(data_dir_files(kk).name, '_');
    dependence = strsplit(data_name{length(data_name)}, 'genes');
    data_name = join(data_name(1:length(data_name)-1), '_');
    data_series{kk-2, 1} = sprintf('%s/%s%s', data_dir_files(kk).folder, data_name{:}, dependence{length(dependence)});
    
    name_arr = strsplit(data_dir_files(kk).name);
    prefix_str = prefix_name;
    for ll=1:length(name_arr),
      if ll==length(name_arr),
        append_name = strsplit(name_arr{ll}, '.');
        if strcmp(append_name(length(append_name)), 'csv'),
          append_name = join(append_name(1:length(append_name)-1), '_');
          append_name = strsplit(append_name{:}, '_');
          append_name = join(append_name(1:length(append_name)-1), '_');
        elseif strcmp(append_name(length(append_name)), 'mat'),
          append_name = {name_arr{ll}};
        else,
          append_name = join(append_name(1:length(append_name)-1), '_');
        end
      else,
        append_name = name_arr{ll};
      end
      prefix_str = sprintf('%s_%s', prefix_str, append_name{:});
    end
    prefix_series{kk-2, 1} = prefix_str;
  end
  data_series = unique(data_series);
  prefix_series = unique(prefix_series);
  medium_series = {};
  medium_series(1:length(prefix_series), 1) = {medium};

end
