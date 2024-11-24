function model=objective_setting_function_v2(model, obj, obj_c, obj_type, model_name),

  %% Testing part
  %% initialize COBRA toolbox
  %addpath('/home/daweilin/cobratoolbox/')
  %run initCobraToolbox()%false)
  %% INPUT files
  %model_path = '/nfs/turbo/umms-csriram/daweilin/data/models/Shen2019.mat'
  %% Load metabolic network model
  %model = load(model_path); % or model_human_duarte.mat
  %fn = fieldnames(model);
  %model = getfield(model, fn{1});
  %model_name = 'Recon1'; %Recon3D%Recon2.2%Recon1
  %% Testing part
  %obj =  {'gh';'atp';'nadh';'nadph';'amet'}
  %obj_c = [0.01, 0.02, 0.01, 0.05, 0.07]
  %obj_type = 'Demand'
  %model_name = 'Recon1'
   
  % getting all compartments 
  %compartments = {};
  %metabolites = {};
  %for i=1:length(model.metabolites),
  %  s = strsplit(model.metabolites{i}, '[');
  %  metabolites{i, 1} = s{1};
  %  ss = strsplit(s{2}, ']');
  %  compartments{i, 1} = ss{1};
  %end
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

  metList = {}
  for iter=1:length(objlist),
      obj = objlist{iter};
      coef = obj_c(iter);
      disp(obj)
      disp(coef)
      % else No objective for optimization
      if strcmp(obj, '')==0,
          if strcmp(obj, 'gh')~=1, % optimize reactions/new demand reactions
              disp('New objective')
              if strcmp('Demand', obj_type),
                  if sum(strfind(obj, '[')),
                      disp('Create a new demand reaction')
                      metList(length(metList)+1, 1) = sprintf('%f %s', coef, obj);
                  else,   
                      disp('Create a new combination demand reaction')
                      obj1 = sprintf('%s%s', obj, compartments{1});
                      obj2 = sprintf('%s%s', obj, compartments{2});
                      obj3 = sprintf('%s%s', obj, compartments{3});
                      obj4 = sprintf('%s%s', obj, compartments{4});
                      obj5 = sprintf('%s%s', obj, compartments{5});
                      obj6 = sprintf('%s%s', obj, compartments{6});
                      obj7 = sprintf('%s%s', obj, compartments{7});
                      inmodel = [findMetIDs(model,obj1)>0, findMetIDs(model,obj2)>0, findMetIDs(model,obj3)>0, findMetIDs(model,obj4)>0, findMetIDs(model,obj5)>0, findMetIDs(model,obj6)>0, findMetIDs(model,obj7)>0];
                      disp('debug!!!!!!!')
                      disp(inmodel)
                      obj_mets = {obj1; obj2; obj3; obj4; obj5; obj6; obj7};
                      obj_mets = obj_mets(find(inmodel));
                      % Add all demand reactions
                      for ii=1:length(obj_mets),
                        obj_mets{ii, 1} = sprintf('%f %s', coef, obj_mets{ii, 1});
                      end
                      join_mets = join(obj_mets, " + ");
                      metList(length(metList)+1, 1) = join_mets;
                  end
              else % optimize biomass optimization
                  disp('Maximize existing reaction')
                  % Change objective reaction if needed
                  exist_ind = find(contains(model.rxns, obj));
                  model.c(exist_ind) = coef;
              end
          else,
              disp('Optimize biomass without using other objectives...')
              % list all the metabolites required for 
              [r_mets, r_coefs] = findMetsFromRxns(model, ori_obj_ind);
              r_mets = r_mets{1}; r_coefs = full(r_coefs{1});
              % Add all demand reactions
              for ii=1:length(r_mets),
                if r_coefs(ii, 1)<0,
                  r_mets{ii, 1} = sprintf('%f %s', -1*coef*r_coefs(ii, 1), r_mets{ii, 1});
                end
              end
              join_mets = join(r_mets, " + ");
              metList(length(metList)+1, 1) = join_mets;
          end
      end
  end
  disp('issue')
  disp(metList)
  % add up the metabolites with corresponding coefficients
  rsym = join(metList, " + ");
  rsym = sprintf('%s -> ', rsym{1});
  rname = sprintf('%s_demand', 'new');
  % add the reaction to the model
  model = addReaction(model, rname, 'reactionFormula', rsym);
  % Change objective reaction if needed
  obj_ind = findRxnIDs(model, rname);
  model.c(obj_ind) = 1;
  disp('Checking objective')
  obj_check = find(model.c~=0);
  model.c(obj_check)
  model.rxns{obj_check}
end
