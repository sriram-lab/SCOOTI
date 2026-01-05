function [fluxstate_uncon, fluxstate_gurobi, grate, geneko_flux, rxnko_growthrate, solverobj, model_out] =  constrain_flux_regulation(model1,onreactions,offreactions,kappa,rho,epsilon,mode,genedelflag,rxndelflag,epsilon2,minfluxflag,CFR_model,extra_weight)


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
    %disp('number of rxns')
    %disp(length(onreactions))
    %disp(length(offreactions))
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

CFR_model_flag = 1;
if (~exist('CFR_model','var')) || (isempty(CFR_model))  

    CFR_model_flag = 0; % by default sum of flux through all ractions is minimized

elseif length(CFR_model)==0,

    CFR_model_flag = 0;
    
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

% get original number of reactions
nrxns = length(model1.rxns);

% stacked CFR models
if CFR_model_flag~=0,
  model = CFR_model;
  extra_weight = extra_weight;
else,
  extra_weight = 1;
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
end
M = 10000;
%epsilon2 = 0;
objpos = find(model.c);
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
    model.obj(colpos) = 1*rho(j)*extra_weight;
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
    model.obj(colpos) = 1*rho(j)*extra_weight;
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
    model.obj(colpos) = -1*kappa(jj)*extra_weight; % minimized
    
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
    model.obj(colpos) = -1*kappa(jj)*extra_weight; % minimized
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
model_out = model;

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


%
close('all')

end

