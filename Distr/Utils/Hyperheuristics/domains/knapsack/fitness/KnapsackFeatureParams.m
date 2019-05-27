function [fP,P,details] = KnapsackFeatureParams(X, fObj, parameters)

% isTransformed = false;
% if ~isempty(varargin), isTransformed = true; M = varargin{1}; W = varargin{2}; end

% printMode = parameters.printMode;
printMode = true;

cutPoint = numel(X)/ 2;
M = X(1:cutPoint);
W = X(cutPoint+1:end);
fObj_ = @(x)fObj(x,M,W);

% Trains selector
if printMode, fprintf('Training over the instances for %d transformation parameters...\t', numel(X)); ...
        fprintf('%.4f\t',X), fprintf('\n'), end
[fP,P,details] = KnapsackTrainSelector(fObj_, parameters);
if printMode, fprintf('Training complete in %d seconds and achieved %.4e of fitness!!!\n', details.time, fP); end

end