function filtered = filter_gene_list(gene_table, model_genes)
  genes = upper(string(gene_table{:,2}));
  filtered = intersect(genes, upper(model_genes));
  filtered = filtered(:);  % Column format
end

