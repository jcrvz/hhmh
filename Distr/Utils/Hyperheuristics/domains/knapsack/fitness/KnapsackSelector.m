function [totalFitness] = KnapsackSelector(X, itemsSet, maxWeightSet, nbRules, nbActions, featureIDs, varargin)    
% Input: X
%
% Input: itemsSet
% It means the wholeData to be used, e.g., {wholeData{trainIndex}}
%
% Input: maxWeightSet
% It means the wholeMaxWeight of the wholeData, e.g., wholeMaxWeight(trainIndex)
%
% Inputs: nbRules, nbActions, featureIDs
% They are the number of rules and actions, and the feature ids to be used
%
% Input: varargin
% M is the first expected varargin and W is the second one ?
% 
% Output: totalFintess
% It is the sum of 


  isTransformed = false;
  if ~isempty(varargin), isTransformed = true; M = varargin{1}; W = varargin{2}; end
  
  % Transforms design vector from the optimizer to Rules+Actions
  [R, A] = getRules(X, nbRules, nbActions); 
  
  if isTransformed
    [solutionTotal] = KnapsackMainSet (itemsSet, maxWeightSet, R, A, featureIDs, M, W);
  else
    [solutionTotal] = KnapsackMainSet (itemsSet, maxWeightSet, R, A, featureIDs);
  end
  totalFitness = sum ([solutionTotal.profit]);
  
end