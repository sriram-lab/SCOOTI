function fixedModel = fixModelRules(model)
    % Build model.rules and model.grRules based on model.rxnGeneMat

    % Extract necessary information
    rxnGeneMat = model.rxnGeneMat; % Binary matrix of reactions and genes
    genes = model.genes; % List of gene IDs in the model

    % Initialize cell arrays for rules and grRules
    numReactions = size(rxnGeneMat, 1);

    fixedModel = model
    fixedModel.rules = cell(numReactions, 1);
    fixedModel.grRules = cell(numReactions, 1);

    % Build the rules and grRules
    for i = 1:numReactions
        genesInReaction = find(rxnGeneMat(i, :));
        if isempty(genesInReaction)
            fixedModel.rules{i} = '';
            fixedModel.grRules{i} = '';
        else
            rule = '';
            grRule = '';
            for j = 1:length(genesInReaction)
                if j > 1
                    rule = [rule ' or '];
                    grRule = [grRule ' or '];
                end
                rule = [rule '(x(' num2str(genesInReaction(j)) '))'];
                grRule = [grRule genes{genesInReaction(j)}];
            end
            fixedModel.rules{i} = rule;
            fixedModel.grRules{i} = grRule;
        end
    end
end

