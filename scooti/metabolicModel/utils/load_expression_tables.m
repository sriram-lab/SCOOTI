function [ups, dws] = load_expression_tables(data_path, upsheet, dwsheet)
  % Loads up- and down-regulated gene tables from CSV, MAT, or XLSX sources.

  if contains(data_path, '.csv')
    % Remove .csv extension
    base_path = erase(data_path, '.csv');
    upfile = sprintf('%s_%s.csv', base_path, upsheet);
    dwfile = sprintf('%s_%s.csv', base_path, dwsheet);
    ups = readtable(upfile);
    dws = readtable(dwfile);

  elseif contains(data_path, '.mat')
    % Load .mat file and use both sheets from the same source
    base_path = erase(data_path, '.mat');
    upfile = sprintf('%s_%s.mat', base_path, upsheet);
    ups_data = load(upfile);
    fn = fieldnames(ups_data);
    ups = ups_data.(fn{1});
    dws = ups;  % assume both are in one struct (adjust if needed)

  else
    % Assume .xlsx or general spreadsheet
    ups = readtable(data_path, 'Sheet', upsheet);
    dws = readtable(data_path, 'Sheet', dwsheet);
  end
end
