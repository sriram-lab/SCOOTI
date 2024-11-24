% Given initial metabolic model and omics input, returns constrained model
% based on INIT

% INPUTS:
% initial_model: unconstrained metabolic model

% alpha: upregulated genes must be greater than mean + alpha * standard
% deviation, downregulated genes must be less than mean - alpha * standard
% deviation

% beta: at least beta percent of genes must meet the criteria described
% just above

% norm_trans: quantile-normalized, log2 normalized omics data

% gene_list: list of genes corresponding to rows in norm_trans and aligning
% with initial_model gene format

% num_threads: number of threads for parallel processing

% OUTPUTS:
% constrained_model: metabolic model constrained by omics gene levels

function [constrained_model] = run_init2(initial_model, upgenes, dwgenes, num_threads, extra_weight)

    %% fix model rule for mapExpressionToReactions
    %initial_model.b = initial_model.b(:, 1);
    %%initial_model = generateRules(initial_model); % didnt work for Recon1 because of the gene ID
    %initial_model = fixModelRules(initial_model);
    %initial_model.c = initial_model.c*0

    % converts gene expression to format for INIT input
    gene_list(1:length(upgenes), 1) = upgenes;
    gene_list(length(upgenes)+1:length(upgenes)+length(dwgenes), 1) = dwgenes;

    gene_expression_input = zeros([length(gene_list) 1]);
    
    gene_expression_input(1:length(upgenes), 1) = 1;
    gene_expression_input(length(upgenes)+1:length(upgenes)+length(dwgenes), 1) = -1;
    
    % multiply by 2 to avoid -1 input, which can cause problems
    gene_expression_input = gene_expression_input .* 1E6;
    
    % prepare data for INIT
    trans_data = struct();
    trans_data.gene = gene_list;
    trans_data.value = gene_expression_input;
    
    trans_rxn_expression = mapExpressionToReactions(initial_model, trans_data);
    
    % sets INIT options
    options = struct();
    options.solver = 'INIT';   
    options.weights = trans_rxn_expression * extra_weight;
    options.timelimit = 45;
    options.printLevel = 1;
    %option.epsilon = 1E-5;
    options.numThreads = num_threads;


    % sets GIMME options
    %options = struct();
    %options.solver = 'GIMME';   
    %options.expressionRxns = trans_rxn_expression * extra_weight;
    %options.threshold = 1E-5;
    %options.timelimit = 45;
    %options.printLevel = 1;
    %options.numThreads = num_threads;
 
    % runs INIT
    constrained_model = INIT(initial_model, options.weights, 1E-8, options.timelimit, 'MILPlog', 1E-5)
    %constrained_model = createTissueSpecificModel(initial_model, options, 1)


end
