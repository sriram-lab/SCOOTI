% Runs CHRR flux sampling, returns average flux vector and all computed fluxes

% NOTE: code based in part on COBRA rMTA tutorial:
% https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorialRMTA.html

% INPUTS
% initial_model: unconstrained metabolic model

% constrained_model: metabolic model constrained by omics

% num_samples: number of flux samples to generate

% OUTPUTS
% flux_avg: vector of average of all sampled fluxes

% flux_all: matrix of all sampled fluxes (total = num_samples)

% rxn_inactive: index of inactive reactions (those not in
% constrained_model)

function [flux_avg, flux_all, rxn_inactive] = run_chrr2(initial_model, constrained_model, num_samples)

    %sampling_options = struct();
    %sampling_options.nPointsReturned = num_samples;   	

    %[~, samples] = sampleCbModel(initial_model, 'sampleFiles', 'CHRR', sampling_options);
        
    sampling_options = struct();
    sampling_options.nPointsReturned = num_samples;   	

    [~, samples] = sampleCbModel(constrained_model, 'sampleFiles', 'CHRR', sampling_options);
        
    % computes sample statistics
    sample_stats = calcSampleStats(samples);

    % converts from reduced index of constrained_model to init_model.rxns index
    idx = zeros(size(samples, 1), 1);

    for i = 1:numel(idx)

        idx(i) = find(strcmp(initial_model.rxns, constrained_model.rxns{i}));
        
    end

    % computes index of inactive reactions
    rxn_inactive = setdiff(1:length(initial_model.rxns), idx); 
  
    % converts sample statistics from constrained model to initial model
    % index
    fields = fieldnames(sample_stats);

    for i = 1:numel(fields)

        temp = sample_stats.(fields{i});
        sample_stats.(fields{i}) = zeros(size(initial_model.rxns));
        sample_stats.(fields{i})(idx) = temp;
        clear temp;

    end
    
    % resizes the flux samples matrix to initial model index 
    temp_2 = samples;
    samples = zeros(size(initial_model.rxns, 1), sampling_options.nPointsReturned);
    samples(idx,:) = temp_2;
    clear temp_2;  

    % outputs all sampled fluxes
    flux_all = samples;

    % computes mean sampled flux
    flux_avg = sample_stats.mean;


end
