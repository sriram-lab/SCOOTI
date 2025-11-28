function model=objective_setting_function(model, obj, obj_c, obj_type, model_name, algorithm),
  % getting all compartments 
  %compartments = {};
  %metabolites = {};
  %for i=1:length(model.metabolites),
  %  s = strsplit(model.metabolites{i}, '[');
  %  metabolites{i, 1} = s{1};
  %  ss = strsplit(s{2}, ']');
  %  compartments{i, 1} = ss{1};
  %end
  % Normalize gene IDs to uppercase without flooding console
  model.genes = upper(model.genes);
  initial_model = model;


  % Create compartment array
  compartments = {'c', 'm', 'n', 'x', 'r', 'g', 'l'};
  for i=1:length(compartments),
    if strcmp(model_name, 'Recon3D'),
      compartments{i} = sprintf('%s%s', '_', compartments{i});
    elseif strcmp(model_name, 'Recon2.2'),
      compartments{i} = sprintf('%s%s', '_', compartments{i});
    else,
      compartments{i} = sprintf('%s%s%s', '[', compartments{i}, ']');
    end
  end
  

  % multi-objective problems
  if isa(obj,'cell')==0,
      objlist = {obj};
      obj_c = [obj_c];
  else,
      objlist = obj;
      obj_c = obj_c;
  end
  % record the original objective
  ori_obj_ind = find(model.c);
  % reset the objective vector
  model.c = model.c*0;
  for iter=1:length(objlist),
      obj = objlist{iter};
      coef = obj_c(iter);
      % else No objective for optimization
      if strcmp(obj, '')==0,
          if strcmp(obj, 'gh')~=1, % optimize reactions/new demand reactions
              %disp('New objective')
              if strcmp('Demand', obj_type),
                  if sum(strfind(obj, '[')),
                      %disp('Create a new demand reaction')
                      demand_rxn = obj;
                      rname = sprintf('%s_demand', obj);
                      rsym = sprintf('%s -> ', obj);
                      model = addReaction(model, rname, 'reactionFormula', rsym);
                      % Change objective reaction if needed
                      obj_ind = find(contains(model.rxns, rname));
                      model.c(obj_ind) = coef;
                  else,   
                      %disp('Create a new combination demand reaction')
                      obj1 = sprintf('%s%s', obj, compartments{1});
                      obj2 = sprintf('%s%s', obj, compartments{2});
                      obj3 = sprintf('%s%s', obj, compartments{3});
                      obj4 = sprintf('%s%s', obj, compartments{4});
                      obj5 = sprintf('%s%s', obj, compartments{5});
                      obj6 = sprintf('%s%s', obj, compartments{6});
                      obj7 = sprintf('%s%s', obj, compartments{7});
                      inmodel = [findMetIDs(model,obj1)>0, findMetIDs(model,obj2)>0, findMetIDs(model,obj3)>0, findMetIDs(model,obj4)>0, findMetIDs(model,obj5)>0, findMetIDs(model,obj6)>0, findMetIDs(model,obj7)>0];
                      %disp('debug!!!!!!!')
                      %disp(inmodel)
                      obj_mets = {obj1; obj2; obj3; obj4; obj5; obj6; obj7};
                      obj_mets = obj_mets(find(inmodel));
                      % Add all demand reactions
                      for ii=1:sum(inmodel),
                          rsym = sprintf('%s -> ', obj_mets{ii});
                          rname = sprintf('%s_demand', obj_mets{ii});
                          model = addReaction(model, rname, 'reactionFormula', rsym);
                          % Change objective reaction if needed
                          obj_ind = findRxnIDs(model, rname);
                          model.c(obj_ind) = coef;
                      end
                  end
              else % optimize biomass optimization
                  %disp('Maximize existing reaction')
                  % Change objective reaction if needed
                  exist_ind = find(contains(model.rxns, obj));
                  model.c(exist_ind) = coef;
              end
          else,
              %disp('Optimize biomass without using other objectives...')

              %% assign objective function
              %if strcmp(model_name, 'Recon3D'),
              %  rxn_name = 'BIOMASS_maintenance';
              %  rf = load('/nfs/turbo/umms-csriram/daweilin/data/BiGG/recon1_biomass_coefficients_for_recon3D.mat')
              %  rf = rf.rf
              %  model = addReaction(model, 'gh_rxn', rf)
              %  gh_ind = find(contains(model.rxns, 'gh_rxn'));
              %  model.c(gh_ind) = coef;
              %elseif strcmp(model_name, 'Recon2.2'),
              %  rxn_name = 'biomass_reaction';
              %  rf = load('/nfs/turbo/umms-csriram/daweilin/data/BiGG/recon1_biomass_coefficients_for_recon2.2.mat')
              %  rf = rf.rf
              %  model = addReaction(model, 'gh_rxn', rf)
              %  gh_ind = find(contains(model.rxns, 'gh_rxn'));
              %  model.c(gh_ind) = coef;
              
              %% original c vector 
              model.c(ori_obj_ind) = coef;

              %% apply recon1 default objective function
              %if strcmp(model_name, 'Recon3D'),
              %  rf = load('/nfs/turbo/umms-csriram/daweilin/data/BiGG/recon1_biomass_coefficients_for_recon3D.mat')
              %  rf = rf.rf
              %  disp(rf)
              %  model = addReaction(model, 'gh_rxn', rf)
              %  gh_ind = find(contains(model.rxns, 'gh_rxn'));
              %  model.c(gh_ind) = coef;
              %elseif strcmp(model_name, 'Recon2.2'),
              %  rf = load('/nfs/turbo/umms-csriram/daweilin/data/BiGG/recon1_biomass_coefficients_for_recon2.2.mat')
              %  rf = rf.rf
              %  model = addReaction(model, 'gh_rxn', rf)
              %  gh_ind = find(contains(model.rxns, 'gh_rxn'));
              %  model.c(gh_ind) = coef;
              %else,
              %  model.c(ori_obj_ind) = coef;
              %end
              %  rxn_name = 'biomass';
              %  gh_ind = find(contains(model.rxns, rxn_name));
              %  model.c(gh_ind) = coef;
              %end
          end
      end
  end
  % Summarize objectives succinctly
  obj_check = find(model.c~=0);
  fprintf('[objective_setting] Objectives set: %d\n', numel(obj_check));
  if ~isempty(obj_check)
    k = min(10, numel(obj_check));
    idx = obj_check(1:k);
    try
      T = table(model.rxns(idx), model.c(idx), 'VariableNames', {'rxn','weight'});
      disp(T)
      if numel(obj_check) > k
        fprintf('[objective_setting] ... and %d more objectives omitted\n', numel(obj_check)-k);
      end
    catch
      % Fallback printing if table fails
      for ii=1:k
        fprintf('  %s : %g\n', model.rxns{idx(ii)}, model.c(idx(ii)));
      end
    end
  end

  %% recover and fix the model and then remove objectives
  %if strcmp(algorithm, 'INIT') | strcmp(algorithm, 'MOOMIN'),
  %  disp('edit initial model')
  %  % fix model rule for mapExpressionToReactions
  %  initial_model.b = model.b(:, 1);
  %  initial_model = fixModelRules(initial_model);
  %  initial_model.c = initial_model.c*0
  %  model = initial_model

  %  unqgenes = unique(model.genes);
  %  model2 = deleteModelGenes(model,unqgenes(1));
  %  disp('Successfully delete a gene 33333')
  %end
  % Ensure objective vector has correct shape and type for downstream Gurobi
  try
    nrxns = length(model.rxns);
    cvec = model.c;
    cvec = full(double(cvec(:))); % force column vector, double
    if length(cvec) < nrxns
      cvec(end+1:nrxns,1) = 0;
    elseif length(cvec) > nrxns
      cvec = cvec(1:nrxns);
    end
    model.c = cvec;
  catch
    % best-effort; leave as-is if any issue
  end
end
