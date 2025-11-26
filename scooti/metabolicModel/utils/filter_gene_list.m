function filtered = filter_gene_list(gene_table, model_genes)
  % Robustly extract a list of genes from a table that may have
  % either one column (named 'upgenes'/'dwgenes') or two+ columns.
  % Falls back to the first available column if the second is absent.

  genes = strings(0,1);
  try
    if istable(gene_table)
      w = width(gene_table);
      if w >= 2
        genes = upper(string(gene_table{:,2}));
      elseif w >= 1
        genes = upper(string(gene_table{:,1}));
      else
        genes = strings(0,1);
      end
    elseif iscell(gene_table)
      genes = upper(string(gene_table(:)));
    else
      genes = upper(string(gene_table));
    end
  catch
    % Fallback: attempt to use the first variable if present
    try
      vars = gene_table.Properties.VariableNames; %#ok<NASGU>
      genes = upper(string(gene_table{:,1}));
    catch
      genes = strings(0,1);
    end
  end

  % Normalize model gene types and intersect
  try
    mg = upper(string(model_genes));
  catch
    mg = upper(string(model_genes(:)));
  end
  filtered = intersect(genes, mg);
  filtered = filtered(:);  % Column format
end
