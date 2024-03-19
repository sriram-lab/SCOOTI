function obj_flux_out=DFAinterface(model_path, obj, obj_type, obj_c, root_path, data_path, sheet_names, upsheet, dwsheet, out_name, ctrl, kappa, geneko_flag, rxnko_flag, medium, media_perturbation, FVAflag)

  sample_name = sprintf('%s%s', upsheet, dwsheet);

  % File name with datetime as prefix
  file_prefix = string(datetime('now','TimeZone','local','Format','MMMdyHHmm'));
  % Output file
  filename = sprintf('%s/[%s]%s', root_path, file_prefix, out_name);
  disp(sprintf('%s %s', 'The DFA result has been saved in', filename));
  % Output file
  excelname = filename;
  
  
  % Load metabolic network model
  model = load(model_path); % or model_human_duarte.mat
  fn = fieldnames(model);
  model = getfield(model, fn{1});
  ms = size(model.b);
  if ms(2)>1,
      model.b = model.b(:,1);
  else,
      model.b = model.b;
  end
  
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % change culture medium KSOM_AA_Jin or DMEMF12_MGSA
  media = readtable('FINAL_MEDIUM_MAP_RECON1.xlsx','Sheet',medium);
  EX_mets = media(:,3);
  EX_rxns = media(:,6);
  [len, ~] = size(EX_rxns);
  m2 = model;
  EX_mets
  for i=1:len,
      ex_rxn_ind = findRxnIDs(model, EX_rxns{i,:});
      if ex_rxn_ind>0;
          m2 = changeRxnBounds(m2, model.rxns{ex_rxn_ind}, EX_mets{i,:}, 'l');
      end
  end
  model = m2;


  % Set objectives for multi- or single-objective problems
  model = objective_setting_function(model, obj, obj_c, obj_type),

  % Save settings in metadata
  metadata.obj = obj;
  metadata.obj_type = obj_type;
  metadata.obj_c = obj_c;
  metadata.output_path = root_path;
  if ctrl==1,
    metadata.input_path = data_path;
  else,
    metadata.input_path = '';
  end
  metadata.file_name = out_name;
  metadata.with_constraint = ctrl;
  metadata.DFA_kappa = kappa;
  metadata.medium = medium;
  metadata.geneko_flag = geneko_flag;
  metadata.rxnko_flag = rxnko_flag;
  metadata.media_perturbation = media_perturbation;
  obj_check = find(model.c~=0); % find objective functions
  if length(obj_check)~=0,
    metadata.objWeights = model.c(obj_check);
    metadata.objRxns = model.rxns{obj_check};
  else
    metadata.objWeights = 0;
    metadata.objRxns = '';
  end
  metadata.model_path = model_path;
  metadata.upStage = upsheet;
  metadata.dwStage = dwsheet;
  
  % convert structure to json files
  encodedJSON = jsonencode(metadata);
  JSONFILE_name= sprintf('%s_metadata.json', excelname);
  fid = fopen(JSONFILE_name,'w');
  fprintf(fid, encodedJSON);
  fclose('all')
  


%% main function for DFA
  function run_simulation(model, excelname),
  
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Initialize progress bar
      progressbarText(0);
  
      % Load essential files as constraints
      T = readtable(data_path,'Sheet',sheet_names);
  
      % Separate values and labels
      % Convert data to num and txt
      num_T = T{:, 2:10};
      txt_T = T{:, 1};
  
      [L, n] = size(num_T);
  
      % Create tables for saving results.
      V_rxn = {}; % reaction deletions
      S_rxn = {};
      V_gene = {}; % gene deletions
      S_gene = {};
  
  
      % RXN DELETIONS
      % Find rxns from the model
      uk = model.rxns;
      disp('The number of reactions:')
      length(uk)
      % Add a row for saving normal GH
      uk{length(uk)+1, 1} = 'WT';
  
      % Create a table for saving results.
      tp1 = "string"; tp2 = "double";
      vn = {}; vn{1} = 'rxns'; vn{2} = sample_name;
      tp = [repelem([tp1, tp2], [1 1])];
      Size = [length(uk) 2];
      V_rxn = table('Size',Size,'VariableTypes',tp,'VariableNames',vn);
      vn = {}; vn{1} = 'mets'; vn{2} = sample_name;
      SizeSlope = [length(txt_T(2:end, 1)) 2];
      S_rxn = table('Size',SizeSlope,'VariableTypes',tp,'VariableNames',vn);
  
      % Save mets name in the first columns
      V_rxn.rxns = uk;
      S_rxn.mets = txt_T(2:end, 1);
  
      disp('Cell type:')
      disp(sheet_names)
      % DFA for embryo
      % kappa = 10
      [~,fluxstate_gurobi,~,~,rc_gurobi,grate_wt,geneko_growthrate,...
      rxnko_growthrate,slope,solverobj,~,~] = flux_activity_coeff2(model,data_path,...
      sheet_names,kappa,1E-3,geneko_flag,rxnko_flag,0,[],[]);
  
      gh_ind = find(contains(model.rxns, 'biomass'));
      gthox_ind = find(contains(model.rxns, '_demand'));
      % model.rxns(gh_ind)
      % res.v(gh_ind)
      z1 = fluxstate_gurobi(gh_ind);
      z2 = sum(fluxstate_gurobi(gthox_ind), 'all');
      % disp('test')
  
      obj_flux_out = [z1 z2]';
  
      % reaction labels
      rxn_labels = model.rxns;
      % additional label for objective value
      rxn_labels{length(model.rxns)+1, 1} = 'Obj';
      
      % Reaction fluxes.
      fluxes = fluxstate_gurobi;
      % Append objective fluxes.
      fluxes(end+1, 1) = grate_wt;

      flux_out = table(model.rxns, fluxstate_gurobi);
      %excelsheet = sprintf('%s_fluxes', sample_name)
      tmp_filename = sprintf('%s_fluxes.csv', excelname)%, excelsheet)
      writetable(flux_out, tmp_filename);
      % writetable(flux_out,excelname,'FileType','spreadsheet','Sheet', excelsheet);
  
      % reduced costs of reactions
      flux_out = table(model.rxns, rc_gurobi);
      excelsheet = sprintf('%s-rc', sample_name)
      tmp_filename = sprintf('%s_%s.csv', excelname, excelsheet)
      writetable(flux_out, tmp_filename);
      % writetable(flux_out,excelname,'FileType','spreadsheet','Sheet', excelsheet);
  
      if rxnko_flag~=0,
          % Find growth rates for gene deletions.
          gh = rxnko_growthrate;
          % Append growth rate for the wild-type strain.
          gh(end+1, 1) = grate_wt;
          % Save result in tb.
          V_rxn.(sample_name) = gh;
          % Save result in tb.
          S_rxn.(sample_name) = slope(:, 1);
  
          % save sheets
          excelsheet = sprintf('%s-rxn', sample_name)
          tmp_filename = sprintf('%s_%s.csv', excelname, excelsheet)
          writetable(V_rxn, tmp_filename);
          % writetable(V_rxn,excelname,'FileType','spreadsheet','Sheet', excelsheet);
          excelsheet = sprintf('%s-rxn-slope', sample_name)
          tmp_filename = sprintf('%s_%s.csv', excelname, excelsheet)
          writetable(S_rxn, tmp_filename);
          % writetable(S_rxn,excelname,'FileType','spreadsheet','Sheet', excelsheet);
      end
  
      % GENE DELETIONS FOR EMBRYO
      % Find rxns from the model
      uk = unique(model.genes);
      uk{length(uk)+1, 1} = 'WT';
  
      % Create a table for saving results.
      tp1 = "string"; tp2 = "double";
      vn = {}; vn{1} = 'genes'; vn{2} = sample_name;
      tp = [repelem([tp1, tp2], [1 1])];
      Size = [length(uk) 2];
      V_gene = table('Size',Size,'VariableTypes',tp,'VariableNames',vn);
  
      vn = {}; vn{1} = 'mets'; vn{2} = sample_name;
      SizeSlope = [length(txt_T(2:end, 1)) 2];
      S_gene = table('Size',SizeSlope,'VariableTypes',tp,'VariableNames',vn);
  
      % Save mets name in the first columns
      V_gene.genes = uk;
      S_gene.mets = txt_T(2:end, 1);
  
      if geneko_flag~=0,
          % Find growth rates for gene deletions.
          gh = geneko_growthrate;
          % Append growth rate for the wild-type strain.
          gh(end+1, 1) = grate_wt;
          % Save result in tb.
          V_gene.(sample_name) = gh;
          % Save result in tb.
          S_gene.(sample_name) = slope(:, 1);
  
          % Save files
          excelsheet = sprintf('%s-gene', sample_name)
          tmp_filename = sprintf('%s_%s.csv', excelname, excelsheet)
          writetable(V_gene, tmp_filename);
          % writetable(V_gene,excelname,'FileType','spreadsheet','Sheet', excelsheet);
          excelsheet = sprintf('%s-gene-slope', sample_name)
          tmp_filename = sprintf('%s_%s.csv', excelname, excelsheet)
          writetable(S_gene, tmp_filename);
          % writetable(S_gene,excelname,'FileType','spreadsheet','Sheet', excelsheet);
      end
      % update progress bar
      progressbarText(i/length(sheet_names));
      
  end
  
  fcn = @run_simulation;
  
  if media_perturbation==1,
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Perturbation of media condition following with the ideas used in egem algorithm
      model2 = model;
      egem = load('./egem/supplementary_software_code.mat');
      mediareactions1 = egem.mediareactions1;
      tmpname = excelname;
      
      delete(gcp('nocreate'))
      parpool(32);
      parfor i = 1:50
          model = model2;
          for kappatype = 1:2
              if kappatype == 1, kappa  = 10; else kappa = 0.01;end
          % delete(gcp('nocreate'))
          % parpool(32);
          % parfor i = 1:50
              excelname = tmpname
              kappa1 = kappa;
              if (kappatype == 2) & (ismember(i,[2,3,5:19])) % trace elements
                  kappa1 = kappa/100;
              elseif (kappatype == 1) & (ismember(i,[1;4])) % glucose or glutamine
                  kappa1 = 3;
              end
              % new suffix of the file names
              if kappatype==1,
                  excelname = sprintf('%s_%s_%s', excelname, 'EXC', mediareactions1{i});
              else,
                  excelname = sprintf('%s_%s_%s', excelname, 'DEP', mediareactions1{i});
              end
              disp(sprintf('%s %s', 'The DFA result has been saved in', excelname));
  
              % change media..
              [ix pos]  = ismember(mediareactions1(i), model.rxns);
              model.lb(pos) = model.lb(pos)*kappa1;
              % run_simulation(model, excelname);
              feval(fcn, model, excelname);
          end
  
          disp(kappatype)
      end
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  else,
      % run_simulation(model, excelname)
      feval(fcn, model, excelname);
  
  end
  % end if media_perturbation







end
