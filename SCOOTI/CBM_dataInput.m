function [data_series, prefix_name, prefix_series, root_path, late_stage, early_stage, simulation, ctrl, medium_series, prefix_pattern_func]=CBM_dataInput(input_data_choose)

  % function to get the name of each column for inferred objective functions
  prefix_pattern_func = @(x) disp(x)

  switch input_data_choose

    %% Significant genes of single cell spatial transcriptome
    case 'CFR_sciSpace'
      % data input setting for single-cell human organ datasets
      data_dir = '/nfs/turbo/umms-csriram/daweilin/data/spatiotemporal/sciSpace/GSE166692_CFR_sampling/sigGenes/Developing-Gut/';
      prefix_name = 'Developing-Gut';
      medium = 'KSOM';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);

      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/sciSpace/Recon1/Developing-Gut/';

    otherwise % CFR without constraint
      data_series = {'/home/daweilin/StemCell/Project_mESC_JinZhang/Jin2021_proteomics/Israel2CBC_proteomics.xlsx'}%;'/home/daweilin/StemCell/Project_mESC_JinZhang/Jin2021_transcriptomics/Gao_2C_BC.xlsx'}
      prefix_name = 'unconstraint'
      prefix_series = {'2CBC_Israel'}%;'2CBC_Gao'}
      %root_path = '/home/daweilin/StemCell/Project_mESC_JinZhang/validation/cancer_model/allmets_models/'
      %root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/pfba/Recon1/KSOM/'
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/pfba/all_mets/DMEMF12/'

      late_stage = 'BC'
      early_stage = '2C'
      simulation = 'model'
      ctrl = 0; % remove constraint
      medium_series = {'DMEMF12'}%;'KSOM'};
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
