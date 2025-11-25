function normalize_Shen2019()
% normalize_Shen2019 — Convert Shen2019.mat string fields to cell arrays of char
% Saves a timestamped backup of Shen2019.mat, then replaces it with normalized version.

  src = fullfile(fileparts(mfilename('fullpath')), 'Shen2019.mat');
  if ~exist(src, 'file')
    error('File not found: %s', src);
  end
  ts = datestr(now, 'yyyymmddHHMMSS');
  bak = [src, '.bak.', ts];
  tmp = [src, '.tmp'];

  S = load(src);
  fns = fieldnames(S);
  model = S.(fns{1});

  % Fields to normalize to 1xN cell array of char
  fields = {'rxns','rxnNames','mets','metNames','genes','grRules','subSystems','metFormulas'};
  for k = 1:numel(fields)
    fld = fields{k};
    if isfield(model, fld)
      v = model.(fld);
      % Convert MATLAB string array to cellstr
      if isstring(v)
        v = cellstr(v);
      end
      % If char array, wrap in cell
      if ischar(v)
        v = {v};
      end
      % If cell, ensure char elements and row shape
      if iscell(v)
        for i = 1:numel(v)
          if isstring(v{i})
            v{i} = char(v{i});
          elseif isnumeric(v{i})
            v{i} = char(string(v{i}));
          elseif iscell(v{i}) && numel(v{i})==1
            w = v{i}{1};
            if ischar(w)
              v{i} = w;
            else
              v{i} = char(string(w));
            end
          elseif ~ischar(v{i})
            try
              v{i} = char(v{i});
            catch
              v{i} = char(string(v{i}));
            end
          end
        end
        v = v(:).'; % force 1xN cell array
      else
        % Last resort: try to convert to cellstr
        try
          v = cellstr(string(v));
          v = v(:).';
        catch
          % leave as is
        end
      end
      model.(fld) = v;
    end
  end

  if isfield(model,'description')
    if isstring(model.description)
      model.description = char(model.description);
    end
  end

  save(tmp, 'model', '-v7');
  movefile(src, bak);
  movefile(tmp, src);
  fprintf('Normalized Shen2019.mat saved. Backup: %s\n', bak);
end

