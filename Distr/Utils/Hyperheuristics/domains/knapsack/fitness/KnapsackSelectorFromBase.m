function [totalFitness] = KnapsackSelectorFromBase(X, itemsSet, maxWeightSet, nbRules, nbActions, featureIDs, baseRules, baseActions, varargin)    
  
  isTransformed = false;
  if ~isempty(varargin), isTransformed = true; M = varargin{1}; W = varargin{2}; end
  
  [R_, A_] = getRules(X, nbRules, nbActions); % Transforms design vector from the optimizer to Rules+Actions
  
  R = [baseRules; R_]; A = [baseActions; A_];
  
  if isTransformed
    [solutionTotal] = KnapsackMainSet (itemsSet, maxWeightSet, R, A, featureIDs, M, W);
  else
    [solutionTotal] = KnapsackMainSet (itemsSet, maxWeightSet, R, A, featureIDs);
  end
  totalFitness = sum ([solutionTotal.profit]);
  
end