% Provides a fitness measure of a a solution
function [fitness] = binpackingFitness(bins)
%  fitness = bins.nbBins;  % Total number of bins used 
%   fitness = mean( bins.freeCapacity ./ bins.capacity );  % Average waste
%   fitness = sum(bins.freeCapacity) / sum(bins.freeCapacity ~= 0);  % Average waste (JCOB)
  fitness = sum(bins.freeCapacity) / length(bins.freeCapacity);  % Average waste (JCOB)
end