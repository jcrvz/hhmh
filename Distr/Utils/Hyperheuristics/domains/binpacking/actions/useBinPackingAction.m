% This function determines the ID of the bin where an item must be packed and packs it
%
% -------------------------------------------------
% Inputs
% -------------------------------------------------
% action:   Integer representing the heuristic to use for selection:
%             1             = Next Fit Heuristic (BP-NF)
%             2             = First Fit Heuristic (BP-FF)
%             3             = Best Fit Heuristic (BP-BF)
%             4             = Worst Fit Heuristic (BP-WF)
%             5             = Almost Worst Fit Heuristic (BP-AWF)
% items_:   Vector containing the data of the items in the instance
% bins:     Structure with the following fields:
%             nbBins        = Current number of bins
%             bin           = Cell array with the information (items) of each bin
%             capacity      = Vector with maximum capacity of each bin
%             freeCapacity  = Vector with remaining capacity of each bin
%
% -------------------------------------------------
% Outputs
% -------------------------------------------------
% items_:   Updated vector containing the data of the items in the instance
% bins:     Updated structure
%
% -------------------------------------------------
% No Example
% -------------------------------------------------

function [items_, bins] = useBinPackingAction(action, items_, bins)
% Select item
[binID] = binpackingStrategy(action, bins, items_(1));
if ~isnan(binID)
    % Check if new bin must be created
    if binID > bins.nbBins
        bins = createNewBin(bins);
    end
    % Assign item to bin and update
    bins.bin{binID} = [bins.bin{binID} items_(1)];
    bins.freeCapacity(binID) = bins.freeCapacity(binID) - items_(1);
end
% Remove item from list
% bins % For debugging purposes. Comment if unrequired
items_ = items_(2:end);
end