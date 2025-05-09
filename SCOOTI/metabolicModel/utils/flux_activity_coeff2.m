function [dynamicmodel,fluxstate_gurobi, fluxstate_gurobi_min, fluxstate_gurobi_max, rc_gurobi, grate_wt, geneko_growthrate, rxnko_growthrate,slope,solverobj,geneko_growthrate_obj,rxnko_growthrate_obj] = flux_activity_coeff2(model, timecourse_metabolomics_datafile,sheetname,kappa,kappa2,genedelflag,rxndelflag,FVAflag,objpos,rxnList)

if (~exist('kappa','var')) || (isempty(kappa))
    kappa = 1;
end

if (~exist('kappa2','var')) || (isempty(kappa2))
    kappa2 = 1E-3;
end

if (~exist('sheetname','var')) || (isempty(sheetname))
    sheetname = 'Sheet1';
end

if (~exist('genedelflag','var')) || (isempty(genedelflag))
    genedelflag = 0;
end

if (~exist('rxndelflag','var')) || (isempty(rxndelflag))
    rxndelflag = 0;
end

if (~exist('FVAflag','var')) || (isempty(FVAflag))
    FVAflag = 0;
end



model_rpmi = model;

if (~exist('objpos','var')) || (isempty(objpos))
    objpos = find(model_rpmi.c);
end

% Specify the default value for rxnList which is a list of all reactions
if (~exist('rxnList','var')) || (isempty(rxnList))
    rxnList = model.rxns;
end

%kappa = 10;
%kappa = 1.5;
%%
% model1 = load('model_MGSA.mat');
% model1 = model1.model_MGSA;
%%
model1 = model_rpmi;
model = model1;
model.A = model1.S;
model.obj = model1.c;
ms = size(model.b);
if ms(2)>1,
    model.rhs = model1.b(:,1);
else,
    model.rhs = model1.b;
end
model.sense =repmat( '=',[size(model1.S,1),1]);
model.lb = model1.lb;
model.ub = model1.ub;
model.vtype = repmat('C',size(model1.S,2),1);
model.modelsense = 'max';
m2 = model;
disp('Successfully get m2')


% load metabolomics data
T = readtable(timecourse_metabolomics_datafile,'Sheet',sheetname);
% Convert data to num and txt
num = T{:, 2:end};
txt = T{:, 1};
manual_matchpos = num(2:end,1:3); % coresponding position in model
timevec =  num(1, 4:end); % time points
%timevec =  num(1, 4:9); % time points
maty = num(2:end,4:end); % metabolomics data
%maty = num(2:end,4:9); % metabolomics data
maty = knnimpute(maty);
slope = zeros(size(maty,1),2);
for i = 1:size(maty,1)
    m1 = mean(maty(i,:));
    [p S] = polyfit(timevec, maty(i,:),1);
    if p(2)==0,
      norm = 1
    else,
      norm = p(2)
    end
    slope(i,:) = p(1)/abs(norm);
    slope(i,:) = slope(i,:)/1000;
%     slope(i,:) = p(1);
end
disp('Successfully get slope')
%[row, col] = find(isnan(slope))
for i = 1:size(maty,1)
    ix0 = (manual_matchpos(:,1) ~= 0);
    
    tol11 = zeros(size(slope(:,1)));
    if ix0(i)
        
        u3pos = manual_matchpos(i,:);
        u3pos(u3pos == 0) = '';
        tol11(i) = slope(i,1); % FLUX ACTIVITY COEFFICIENT
        epsilon2 = tol11(i);
        
        rowpos = u3pos;
        colpos = size(m2.A,2) + 1;
        colpos1 = colpos + 1;
        m2.rhs(rowpos) = -epsilon2;
        
        
        % metrow = tol - alpha + beta
        
        m2.A(rowpos,colpos) = 1;
        m2.A(rowpos,colpos1) = -1;
        m2.vtype(colpos) = 'C';
        
        % set si to be positive
        m2.lb(colpos1) = 0;
        m2.ub(colpos1) = 1000;
        m2.lb(colpos) = 0;
        m2.ub(colpos) = 1000;
        m2.obj(colpos1) = -1*kappa; % minimized
        m2.obj(colpos) = -1*kappa; % minimized               
    end
end
disp('Successfully process the first for-loop')
epsilon = 0;
%kappa2 = 5.5E-3;
for jj = 1:length(model.rxns)
    if model.c(jj)==0
        rxnpos = jj;
        %        xi + si >= -eps2
        %     si >= 0
        %     rho(ri + si)
        % constraint 1
        rowpos = size(m2.A,1) + 1;
        colpos = size(m2.A,2) + 1;
        disp(size(m2.rhs))
        m2.A(rowpos,rxnpos) = 1;
        m2.A(rowpos,colpos) = 1;
        m2.rhs(rowpos) = -epsilon;
        m2.sense(rowpos) = '>';
        % set si to be positive
        m2.lb(colpos) = 0;
        m2.ub(colpos) = 1000;
        m2.obj(colpos) = -1*kappa2; % minimized
    
        % constraint 2
        %     xi - ri <= eps2
        %     ri >= 0
        % new row and column
        rowpos = size(m2.A,1) + 1;
        colpos = size(m2.A,2) + 1;
        m2.A(rowpos,rxnpos) = 1;
        m2.A(rowpos,colpos) = -1;
        m2.rhs(rowpos) = epsilon;
        m2.sense(rowpos) = '<';
        % set ri to be positive
        m2.lb(colpos) = 0;
        m2.ub(colpos) = 1000;
        m2.obj(colpos) = -1*kappa2; % minimized
    end
end
disp('objective setting')
% find(model.c==1)
% find(m2.obj==1)
disp('Successfully process the second for-loop')
m2.A = sparse(m2.A);
m2.vtype = repmat('C',size(m2.A,2),1);
params.outputflag = 0;
%[row, col] = find(isnan(m2.rhs))
solg1 = gurobi(m2,params);
grate_wt =  sum(solg1.x(objpos));
solverobj = solg1;
dynamicmodel = m2;

% report reduced costs for reactions without perturbations
rc_gurobi = solg1.rc(1:length(model.rxns));
% report fluxes without perturbations
try
    fluxstate_gurobi = solg1.x(1:length(model.rxns));
catch
    fluxstate_gurobi = NaN;
end

% FVA with gurobi
if FVAflag,
    disp('Processing FVA...')
    m3 = m2;
    m3.lb(find(m3.c)) = solg1.x(find(m3.c));
    rxn_ids = find(strcmp(rxnList, m3.rxns));
    for ii=1:length(rxn_ids),
        disp(ii)
        rr = rxn_ids(ii);
        m3.c = m3.c*0;
        m3.c(rr,1) = 1;
        m3.modelsense = 'min';
        solg_min = gurobi(m3,params);
        m3.modelsense = 'max';
        solg_max = gurobi(m3,params);
        fluxstate_gurobi_min(ii,1) = solg_min.x(ii);
        fluxstate_gurobi_max(ii,1) = solg_max.x(ii);
    end
else
    fluxstate_gurobi_min(1:length(m2.rxns),1) = zeros(length(m2.rxns), 1);
    fluxstate_gurobi_max(1:length(m2.rxns),1) = zeros(length(m2.rxns), 1);
end


% optional. do gene deletion analysis
if genedelflag
    unqgenes = unique(m2.genes);
    for kk = 1:length(unqgenes),
        model2 = deleteModelGenes(m2,unqgenes(kk));
        solg1 = gurobi(model2,params);
        geneko_growthrate(kk,1) = sum(solg1.x(find(model_rpmi.c)));
        geneko_growthrate_obj(kk,1) = solg1.objval;
        
        disp(kk)
        
    end
else
    geneko_growthrate = [];
    geneko_growthrate_obj = [];
end

if rxndelflag
    
    for kk = 1:length(model_rpmi.rxns),
        model2 = m2; model2.lb(kk) = 0; model2.ub(kk) = 0;
        solg1 = gurobi(model2,params);
        rxnko_growthrate(kk,1) = sum(solg1.x(find(model_rpmi.c)));
        rxnko_growthrate_obj(kk,1) = solg1.objval;
        
        disp(kk)
    end
    
else
    rxnko_growthrate = 0;
    rxnko_growthrate_obj = [];
end

close('all')

end

    
