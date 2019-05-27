% X = Design vector used by optimizer
% itemsSet = Set of instances
% defaultWeightSet = Set of default weight for each container in the set of instances
% nbRules = Number of rules in the selector
% nbActions = Number of possible heuristics in the selector

function [totalFitness] = BinPackingSelector(X, itemsSet, defaultWeightSet, nbRules, nbActions, featureIDs, varargin)
isTransformed = false;
isPCA = false;

if ~isempty(varargin)
    if length(varargin) == 2
        isTransformed = true; M = varargin{1}; W = varargin{2};
    elseif length(varargin) == 3
        isPCA = true;
        mu = varargin{1}; sigma = varargin{2}; U = varargin{3};
    end
end

[R, A] = getRules(X, nbRules, nbActions); % Transforms design vector from the optimizer to Rules+Actions

if isTransformed
    [solutionTotal] = BinPackingMainSet (itemsSet, defaultWeightSet, R, A, featureIDs, M, W);
elseif isPCA
    [solutionTotal] = BinPackingMainSet (itemsSet, defaultWeightSet, R, A, featureIDs, mu, sigma, U);
else
    [solutionTotal] = BinPackingMainSet (itemsSet, defaultWeightSet, R, A, featureIDs);
end

totalFitness = sum (solutionTotal.fitness);
end