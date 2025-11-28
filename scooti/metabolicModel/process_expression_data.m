function [ups_table, dws_table] = process_expression_data(ups, dws, model_name)
  % Choose BiGG mapping path based on model
  switch model_name
    case 'Recon3D'
      map_path = './SCOOTI/SCOOTI/metabolicModel/GEMs/Recon3D_genes.json';
    case 'Recon2.2'
      map_path = './SCOOTI/SCOOTI/metabolicModel/GEMs/Recon2.2_symbol_to_hgnc.json';
    case 'Recon1'
      % Optional mapping for Recon1 if available; otherwise fall back to uppercase only
      cand = './SCOOTI/SCOOTI/metabolicModel/GEMs/Recon1_symbol_to_hgnc.json';
      if exist(cand, 'file')
        map_path = cand;
      else
        warning('Recon1 mapping JSON not found; using uppercase gene symbols without conversion.');
        ups_table = ups; dws_table = dws; return;
      end
    otherwise
      warning('Model name not recognized for BiGG mapping. Skipping annotation.');
      ups_table = ups; dws_table = dws; return;
  end

  % Ensure gene names are uppercase
  ups{:,2} = upper(ups{:,2});
  dws{:,2} = upper(dws{:,2});

  % Convert using mapping
  up_genes = BiGG_to_annotation(ups{:,2}, map_path, 'to_BiGG');
  dw_genes = BiGG_to_annotation(dws{:,2}, map_path, 'to_BiGG');

  % Wrap back into tables
  ups_table = array2table([(1:length(up_genes))', up_genes], ...
      'VariableNames', {'Index', 'Gene'});
  dws_table = array2table([(1:length(dw_genes))', dw_genes], ...
      'VariableNames', {'Index', 'Gene'});

end
