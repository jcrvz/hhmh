% This function determines the ID of the bin where an item must be packed
%
% -------------------------------------------------
% Inputs
% -------------------------------------------------
% mode:     Integer representing the heuristic to use for selection:
%             1             = Next Fit Heuristic (BP-NF)
%             2             = First Fit Heuristic (BP-FF)
%             3             = Best Fit Heuristic (BP-BF)
%             4             = Worst Fit Heuristic (BP-WF)
%             5             = Almost Worst Fit Heuristic (BP-AWF)
% bins:     Structure with the following fields:
%             nbBins        = Current number of bins
%             bin           = Cell array with the information (items) of each bin
%             capacity      = Vector with maximum capacity of each bin
%             freeCapacity  = Vector with remaining capacity of each bin
% item:     Data (e.g. weight) of the item that will be packed
%
% -------------------------------------------------
% Outputs
% -------------------------------------------------
% binID:    ID of the bin where the next item will be packed. If higher than
%           nbBins, it means that a new bin should be opened
%
% -------------------------------------------------
% Example
% -------------------------------------------------
% ti = [4, 8, 5, 1, 7, 6, 1, 4, 2, 2]; tb = struct('nbBins',1,'freeCapacity', [10]);
% binpackingStrategy (5,tb,ti(1))
% returns:
%   ans =  1
function [binID] = binpackingStrategy(mode, bins, item)
% Begin item validation
if all(item > bins.capacity)
    % Item is invalid
    binID = NaN;
    warning("Invalid item detected (bigger than all bins)... Discarding the item ...")
    return
end
% End item validation
switch mode
    case 1      % Next fit
        if item <= bins.freeCapacity(bins.nbBins) % Check current bin
            binID = bins.nbBins;
        else
            binID = bins.nbBins + 1; % No bin was found so assign a new one
        end
    case 2     % First fit
        for idx = 1 : bins.nbBins
            if item <= bins.freeCapacity(idx)
                binID = idx;
                return; %  Appropriate bin was found. Exit loop and return data
            end
        end
        binID = idx + 1; % No bin was found so assign a new one
    case 3     % Best fit
        feasibleBins = bins.freeCapacity >= item;
        if any(feasibleBins)
            [~, binID] = min( bins.freeCapacity(feasibleBins) ); % Gives biggest free capacity
            binIndex = 1:bins.nbBins;
            validBins = binIndex(feasibleBins);
            binID = validBins(binID);
        else
            binID = bins.nbBins + 1; % No bin was found so assign a new one
        end
    case 4     % Worst fit
        [binRemainingCapacity, binID] = max( bins.freeCapacity ); % Gives biggest free capacity
        if item > binRemainingCapacity % Check if item does not fit
            binID = bins.nbBins + 1; % No bin was found so assign a new one
        end
    case 5     % Almost-Worst Fit
        feasibleBins = bins.freeCapacity >= item;
        nbFeasibleBins = sum(feasibleBins); % Checks how many are True
        binIndex = 1 : bins.nbBins;
        if nbFeasibleBins == 1 % Only one valid bin
            binID = binIndex(feasibleBins); % Indexing returns a scalar
        elseif nbFeasibleBins > 1 % Two or more available
            [sortedCapacity, sortedIndex] = sort( bins.freeCapacity(feasibleBins), 'descend' );
            validBins = binIndex(feasibleBins);
            binID = validBins(sortedIndex(2));
        else
            binID = bins.nbBins + 1; % No bin was found so assign a new one
        end
    otherwise
        binID = NaN;
        error("Invalid heuristic requested (out of range). Aborting! ...")
end
end