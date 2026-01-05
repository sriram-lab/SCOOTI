function [data_series, prefix_series, medium_series] = split_series(data_series_str, prefix_series_str, medium_series_str)
  data_series = strsplit(data_series_str, ',');
  prefix_series = strsplit(prefix_series_str, ',');
  medium_series = strsplit(medium_series_str, ',');
end
