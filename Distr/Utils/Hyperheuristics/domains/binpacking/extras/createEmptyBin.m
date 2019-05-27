function [bin] = createEmptyBin(defaultWeight)
bin = struct('nbBins',0, 'defaultCapacity', defaultWeight); % Creates empty solution
bin = createNewBin (bin); % Adds the first container
end