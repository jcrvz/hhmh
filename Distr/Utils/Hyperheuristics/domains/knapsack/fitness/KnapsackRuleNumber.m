function [fP,P,details] = KnapsackRuleNumber(X, fObj, parameters)

% isTransformed = false;
% if ~isempty(varargin), isTransformed = true; M = varargin{1}; W = varargin{2}; end

% printMode = parameters.printMode;
printMode = true;

nbRules = round(X);
parameters.numRules = nbRules; % 

fObj_ = @(x)fObj(x,nbRules);

% Trains selector
if printMode, fprintf('Training over the instances for %d rules...\n', nbRules); end
[fP,P,details] = KnapsackTrainSelector(fObj_, parameters);
if printMode, fprintf('Training complete in %d seconds and achieved %.4e of fitness!!!\n', details.time, fP); end

end