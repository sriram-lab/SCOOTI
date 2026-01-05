function geneConverted = BiGG_to_annotation(geneList, gene_info_path, mode),
  
  % Converted BiGG gene ID to gene annotations
  %
  % Arguments
  % ---------
  %   geneList (cell array): BiGG gene IDs or gene names
  %   gene_info_path (string): path to access the information of genes
  %   mode (string): 'to_BiGG' or 'to_annot'
  
  % read files
  fileName = gene_info_path; % filename in JSON extension
  fid = fopen(fileName); % Opening the file
  raw = fread(fid,inf); % Reading the contents
  str = char(raw'); % Transformation
  fclose(fid); % Closing the file
  data = jsondecode(str); % Using the jsondecode function to parse JSON from string

  % get basic info (support both Recon3D-style and simple arrays regardless of filename)
  if isstruct(data) && isfield(data, 'results')
    try
      gene_names = {data.results.name};
      gene_BiGGID = {data.results.bigg_id};
    catch
      error('Unexpected JSON schema in results for %s', fileName);
    end
  else
    if isfield(data, 'symbol') && isfield(data, 'hgnc')
      gene_names = data.symbol;
      gene_BiGGID = data.hgnc;
    else
      error('Unsupported JSON schema in %s (expected results[] or symbol/hgnc arrays).', fileName);
    end
  end

  % Normalize to strings and uppercase for robust matching (return original strings)
  gene_names = string(gene_names(:));
  gene_BiGGID = string(gene_BiGGID(:));
  names_upper = upper(strtrim(gene_names));
  ids_upper   = upper(strtrim(gene_BiGGID));

  % Initiating conversion...
  geneAnnotations = [];
  geneBiGGs = [];
  geneConverted = [];
  sg = size(geneList);
  sg = sg(1);
  for i=1:sg,
    q = upper(string(geneList(i,1)));
    if strcmp(mode, 'to_annot'),
      % map BiGG ID -> symbol
      ind = find(strcmp(q, ids_upper));
      if ~isempty(ind)
        geneConverted = [geneConverted, cellstr(gene_names(ind))];
      end
    else
      % map symbol -> BiGG ID
      ind = find(strcmp(q, names_upper));
      if ~isempty(ind)
        geneConverted = [geneConverted, cellstr(gene_BiGGID(ind))];
      end
    end
  end
  % transpose
  geneConverted = transpose(geneConverted);
end


