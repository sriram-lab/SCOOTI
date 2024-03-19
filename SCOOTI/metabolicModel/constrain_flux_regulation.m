function [fluxstate_uncon, fluxstate_gurobi, grate, geneko_flux, rxnko_growthrate, solverobj, grate_wt_min, grate_wt_max] =  constrain_flux_regulation(model1,onreactions,offreactions,kappa,rho,epsilon,mode,genedelflag,rxndelflag,epsilon2,minfluxflag,FVAflag)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist('mode','var')) || (isempty(mode))
    mode = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('on_off_reactions:')
%disp(onreactions)
%disp(offreactions)
if mode == 0 % genes
    [~,~,onreactions,~] =  deleteModelGenes(model1, onreactions);
    [~,~,offreactions,~] =  deleteModelGenes(model1, offreactions);
    disp('number of rxns')
    disp(length(onreactions))
    disp(length(offreactions))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist('epsilon','var')) || (isempty(epsilon))
    
    epsilon = ones(size(onreactions))*1E-3;
    
end

if numel(epsilon) == 1
    
    epsilon = repmat(epsilon, size(onreactions));
    
end


if (~exist('rho','var')) || (isempty(rho))
    
    rho = repmat(1, size(onreactions));
    
end

if numel(rho) == 1
    
    rho  = repmat(rho, size(onreactions));
    
end

if (~exist('kappa','var')) || (isempty(kappa))
    
    kappa = repmat(1, size(offreactions));

end

if numel(kappa) == 1
    
    kappa  = repmat(kappa, size(offreactions));
    
end


if (~exist('epsilon2','var')) || (isempty(epsilon2))
    
    epsilon2 = zeros(size(offreactions));
end

if (~exist('minfluxflag','var')) || (isempty(minfluxflag))
    
    minfluxflag = true; % by default sum of flux through all ractions is minimized
    
end

kappa1 = 1E-4; % minfluxflag
params.outputflag = 0;
% params.iterationlimit = 1E9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if minfluxflag
kappa = [kappa(:); ones(size(setdiff(model1.rxns, offreactions)))*kappa1]; % minimize flux through all the reactions. PFBA
epsilon2 = [epsilon2; zeros(size(setdiff(model1.rxns, offreactions)))];
offreactions = [offreactions(:); setdiff(model1.rxns, offreactions)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert to gurobi format.
% model1 = load('model_MGSA.mat');
% model1 = model1.model_MGSA;
model = model1;
model.A = model1.S;
model.obj = model1.c;
model.rhs = model1.b;
if exist('model1.csense','var') && ~isempty(model1.csense)
    model.sense = model1.csense;
    model.sense(ismember(model.sense,'E')) = '=';
    model.sense(ismember(model.sense,'L')) = '<';
    model.sense(ismember(model.sense,'G')) = '>';
else
    model.sense =repmat( '=',[size(model1.S,1),1]);
end
model.lb = model1.lb;
model.ub = model1.ub;
model.vtype = repmat('C',size(model1.S,2),1);
model.modelsense = 'max';
nrows = size(model.A,1);
ncols = size(model.A,2);
M = 10000;
%epsilon2 = 0;
objpos = find(model.c);
nrxns = length(model.rxns);
% solve fluxes without any constraints
solg_uncon = gurobi(model,params);
fluxstate_uncon = solg_uncon.x(1:nrxns);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maximize the number of reactions with proteomic or transcriptomic evidence that are ON/up-regulated.
for j = 1:length(onreactions)
    
    rxnpos = find(ismember(model1.rxns,onreactions(j)));
    
    %         xi - (eps + M)ti >= -M
    % ti = 0 or 1.
    
    rowpos = size(model.A,1) + 1;
    colpos = size(model.A,2) + 1;
    
    model.A(rowpos,rxnpos) = 1;
    model.A(rowpos,colpos) = -(1*epsilon(j) + M);
    model.rhs(rowpos) = -M;
    model.sense(rowpos) = '>';
    model.vtype(colpos) = 'B';
    model.obj(colpos) = 1*rho(j);
    model.lb(colpos) = 0;
    model.ub(colpos) = 1;
    
    % xi + (eps + M)ri <= M
    % ri = 0 or 1.
    
    rowpos = size(model.A,1) + 1;
    colpos = size(model.A,2) + 1;
    
    model.A(rowpos,rxnpos) = 1;
    model.A(rowpos,colpos) = (1*epsilon(j) + M);
    model.rhs(rowpos) = M;
    model.sense(rowpos) = '<';
    model.vtype(colpos) = 'B';
    model.obj(colpos) = 1*rho(j);
    model.lb(colpos) = 0;
    model.ub(colpos) = 1;
    
end

%  constraints for off reactions. their flux is minimized.
% soft constraints - can be violated if neccesary - i.e some reactions
% can have some flux.  higher magnitude higher penalty

for jj = 1:length(offreactions)
    rxnpos = find(ismember(model1.rxns,offreactions(jj)));
    %        xi + si >= -eps2
    %     si >= 0
    %     rho(ri + si)
    % constraint 1
    rowpos = size(model.A,1) + 1;
    colpos = size(model.A,2) + 1;
    model.A(rowpos,rxnpos) = 1;
    model.A(rowpos,colpos) = 1;
    model.rhs(rowpos) = -epsilon2(jj);
    model.sense(rowpos) = '>';
    model.vtype(colpos) = 'C';
    % set si to be positive
    model.lb(colpos) = 0;
    model.ub(colpos) = 1000;
    model.obj(colpos) = -1*kappa(jj); % minimized
    
    % constraint 2
    %     xi - ri <= eps2
    %     ri >= 0
    % new row and column
    rowpos = size(model.A,1) + 1;
    colpos = size(model.A,2) + 1;
    model.A(rowpos,rxnpos) = 1;
    model.A(rowpos,colpos) = -1;
    model.rhs(rowpos) = epsilon2(jj);
    model.sense(rowpos) = '<';
    model.vtype(colpos) = 'C';
    % set ri to be positive
    model.lb(colpos) = 0;
    model.ub(colpos) = 1000;
    model.obj(colpos) = -1*kappa(jj); % minimized
end

solg1 = gurobi(model,params);
try
    fluxstate_gurobi = solg1.x(1:nrxns);
    grate = sum(solg1.x(objpos));
    solverobj = solg1.objval;
catch
    fluxstate_gurobi = NaN;
    grate = NaN;
    solverobj = NaN;
end


% optional. do gene deletion analysis
if genedelflag
    geneko_flux(:,1) = solg1.x(1:nrxns);
    unqgenes = unique(model.genes);
    for kk = 1:length(unqgenes),
        model2 = deleteModelGenes(model,unqgenes(kk));
        solg1 = gurobi(model2,params);
        geneko_growthrate(kk,1) = sum(solg1.x(find(model1.c)));
        geneko_growthrate_obj(kk,1) = solg1.objval;
        % record fluxes after knockouts
        disp(solg1.objval)
        fluxstate_ko = solg1.x(1:nrxns);
        geneko_flux(:,kk+1) = fluxstate_ko;
        disp(kk)
    end
else
    geneko_growthrate = [];
    geneko_growthrate_obj = [];
    geneko_flux = [];
end

if rxndelflag
    
    for kk = 1:length(model.rxns),
        model2 = model; model2.lb(kk) = 0; model2.ub(kk) = 0;
        solg1 = gurobi(model2,params);
        rxnko_growthrate(kk,1) = sum(solg1.x(find(model1.c)));
        rxnko_growthrate_obj(kk,1) = solg1.objval;
        
        disp(kk)
    end
    
else
    rxnko_growthrate = 0;
    rxnko_growthrate_obj = [];
end


% optional. do FVA to get ranges of variables that achieve the same objective values
%disp('Processing FVA...')

grate_wt_min = 0;
grate_wt_max = 0;
if FVAflag,

  % FVA with gurobi
  disp('Processing FVA...')

  %params.outputflag = 0;
  %model = load('cancer_model.mat')
  %model = model.model;
  %objind = find(contains(model.rxns, 'biomass'));
  %model.c(objind) = 1;
  %model1 = model;
  %model.A = model1.S;
  %model.obj = model1.c;
  %model.rhs = model1.b;
  %model.lb = model1.lb;
  %model.ub = model1.ub;


  objind = find(model.obj~=0);
  %model.lb(objind) = solg1.x(objind);
  %model.ub(objind) = solg1.x(objind);
  %model.sense =repmat( '=',[size(model1.S,1),1]);
  %model.vtype = repmat('C',size(model1.S,2),1);
  %model.modelsense = 'min';
  %solg2 = gurobi(model, params);
  
  m3 = model;

  %m3.rhs(find(model.obj~=0)) = solg1.x(find(model.obj~=0));
  %m3.sense(find(model.obj~=0)) = repmat( '=',[length(find(model.obj~=0)),1]);

  m3.lb(find(model.obj~=0)) = solg1.x(find(model.obj~=0));
  m3.ub(find(model.obj~=0)) = solg1.x(find(model.obj~=0));
  grate_wt = 0;
  grate_wt_min = 0;
  grate_wt_max = 0;
  infct = 0;
  for kk = 1:nrxns,%length(model.rxns),
      m3.obj = m3.obj*0;
      m3.obj(kk,1) = 1;
      m3.modelsense = 'min';
      solg_min = gurobi(m3,params);
      disp(solg_min)
      %disp(solg_min.x(kk))
      %disp(solg_min.objval)
      disp('get min')

      m3.obj(kk,1) = 1;
      m3.modelsense = 'max';
      solg_max = gurobi(m3,params);
      disp(solg_max)
      try
        min_v = solg_min.x(kk);
      catch
      %if strcmp(solg_min.status, 'INFEASIBLE'),
      %  infct = infct+1;
        min_v = solg1.x(kk);
      end
      try
        max_v = solg_max.x(kk);
      catch
        max_v = solg1.x(kk);
      %else,
      %  disp(solg_min)
      %  %disp(solg_min.x(kk))
      %  %disp(solg_min.objval)
      %  disp('get min')
      %  m3.modelsense = 'max';
      %  solg_max = gurobi(m3,params);
      %  min_v = solg_min.x(kk);
      %  max_v = solg_max.x(kk);
      end
      grate_wt_min(kk,1) = min_v;%min(min_v, max_v);
      grate_wt_max(kk,1) = max_v;%max(min_v, max_v);
      %grate_wt(kk,1) = solg1.x(kk);
  end
  infct
end
%
close('all')

end

