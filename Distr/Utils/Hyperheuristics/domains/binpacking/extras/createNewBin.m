% This function adds a new bin to the structure containing all bins 
% 
% -------------------------------------------------
% Inputs 
% -------------------------------------------------
% bins:     Structure with at least the following fields: 
%             nbBins          = Current number of bins
%             bin             = Cell array with the information (items) of each bin 
%             capacity        = Vector with maximum capacity of each bin
%             freeCapacity    = Vector with remaining capacity of each bin
%             defaultCapacity = Default value for all bins (assuming a single capacity for all)
%
% -------------------------------------------------
% Outputs
% -------------------------------------------------
% bins:     Modified structure with the new bin
%
% -------------------------------------------------
% Example 
% -------------------------------------------------
% tb = struct('nbBins',1,'freeCapacity', [10],'bin',{{[]}}, 'defaultCapacity', 10, 'capacity', 10);
% tb = createNewBin (tb)
% returns:
%  scalar structure containing the fields:
%  
%      nbBins =  2
%      freeCapacity =
%          
%            10   10                 
%      bin = 
%      {
%            [1,1] = [](0x0)
%            [1,2] = [](0x0)
%      }      
%      defaultCapacity =  10
%      capacity = 
%      
%             10   10


function [bins] = createNewBin(bins)
  newBinID = bins.nbBins + 1;
  bins.bin{newBinID} = [];
  bins.capacity(newBinID) = bins.defaultCapacity;
  bins.freeCapacity(newBinID) = bins.capacity(newBinID);
  bins.nbBins = newBinID;
end