function [data_series, prefix_name, prefix_series, root_path, late_stage, early_stage, simulation, ctrl, medium_series, prefix_pattern_func]=CBM_dataInput(input_data_choose)

  % function to get the name of each column for inferred objective functions
  prefix_pattern_func = @(x) disp(x)

  switch input_data_choose

%% Flux prediction of gene knockout strains based on bulk transcriptomics (technically the same as CFR_bulkomics)
    case 'CFR_knockout'
      % data input setting for single-cell human organ datasets

      data_dir = './example_omics_data/CFR_knockout/';
      prefix_name = 'ESC';
      medium = 'DMEMF12';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);

      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = './example_fluxPrediction/CFR_knockout/';
   

%% Flux prediction based on bulk transcriptomics/proteomics
    case 'CFR_bulkOmics'
      data_dir = './example_omics_data/CFR_bulkOmics/';
      prefix_name = 'prolif-qui';
      medium = 'DMEMF12';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);

      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = './fluxPrediction/CFR_bulkOmics/';
   
%% Flux prediction based on single-cell datasets
    case 'CFR_scOmics'
      % input data
      data_dir = './example_omics_data/CFR_scOmics/';
      prefix_name = 'scCellCycle';
      medium = 'DMEMF12';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);
      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = './example_fluxPrediction/CFR_scOmics/';
      % function to get cell-specific objectives
      prefix_pattern_func = @prefix_pattern_scCellCycleHela

    
%% Flux prediction based on metabolomics datasets
    case 'DFA_bulkOmics'
      data_series = {'./example_omics_data/Jin_Time_Course_Metabolomics.xlsx';'./example_omics_data/Jin_Metabolomics_DFA.xlsx'};
      prefix_name = '1C2C';
      prefix_series = {'E1C-L2C';'ES-ES'};
      root_path = './example_fluxPrediction/DFA_bulkOmics/';
      late_stage = '2C';
      early_stage = '1C';
      simulation = 'DFA';
      ctrl = 1; %set to 1 if we want to apply constraints
      medium_series = {'KSOM';'DMEMF12'};

%% Flux prediction with pFBA without constraint
    otherwise 
      data_series = {'./example_omics_data/Israel2CBC_proteomics.xlsx'};
      prefix_name = 'unconstraint';
      prefix_series = {'2CBC_Israel'};
      root_path = './example_fluxPrediction/unconstraint_models/';
      late_stage = 'BC'
      early_stage = '2C'
      simulation = 'model'
      ctrl = 0; % remove constraint
      medium_series = {'DMEMF12';'KSOM'};
  end
end



function prefix_pattern_name=prefix_pattern_scEmbryo(prefix)
  prefix_pattern_name = strsplit(prefix, '_');
  prefix_pattern_name = join(prefix_pattern_name(2:end), '_');
end


function prefix_pattern_name=prefix_pattern_scCellCycleHela(prefix)
  prefix_pattern_name = strsplit(prefix, '_');
  prefix_pattern_name = join(prefix_pattern_name(4:end), '_');
end
