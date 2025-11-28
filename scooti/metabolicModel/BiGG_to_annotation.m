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
      % Fallback if results layout differs
      error('Unexpected JSON schema in results for %s', fileName);
    end
  else
    % Expect parallel arrays: symbol (names), hgnc (BiGG IDs)
    if isfield(data, 'symbol') && isfield(data, 'hgnc')
      gene_names = data.symbol;
      gene_BiGGID = data.hgnc;
    else
      error('Unsupported JSON schema in %s (expected results[] or symbol/hgnc arrays).', fileName);
    end
  end

  % Initiating conversion...
  geneAnnotations = [];
  geneBiGGs = [];
  geneConverted = [];
  sg = size(geneList);
  sg = sg(1);
  for i=1:sg,
    if strcmp(mode, 'to_annot'),
      ind = find(strcmp(geneList(i,1), gene_BiGGID));
      if length(ind)>0,
        geneConverted = [geneConverted, gene_names(ind)];
      end
    else,
      ind = find(strcmp(geneList(i,1), gene_names));
      if length(ind)>0,
        geneConverted = [geneConverted, gene_BiGGID(ind)];
      end
    end
  end
  % transpose
  geneConverted = transpose(geneConverted);
end



