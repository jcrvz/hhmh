% Function for the Knapsack (KP) problem (solves a set of instances)
% -------------------------------------------------
% Inputs 
% -------------------------------------------------
% itemsSet:     Cell array. Each position contains a data matrix of items, where: Row: Item. Columns: Weight, Profit
% maxWeightSet: Vector with maximum weight for each instance 
% rules:        Matrix containing the set of rules. One row per rule and one column per feature
% actions:      Column vector with action indices, starting at one
%
% -------------------------------------------------
% Outputs 
% -------------------------------------------------
% solutionTotal: Structure containing the solution. It has the following fields: fitness, weight, items. 
%                The first two (arrays) give the total fitness and total weight for each instance. The 
%                last one gives the packed items for each instance. 
%
% -------------------------------------------------
% Example 
% -------------------------------------------------
% [wholeData, wholeMaxWeight, sumItems] = loadInstances(0);
% mw = 20;
% ts = [0.5*ones(1,7)];
% ta = [1];
% KnapsackMainSet(wholeData, wholeMaxWeight, ts, ta);
% 
% Returns:
%   scalar structure containing the fields:
%
%    fitness =
%
%     Columns 1 through 16:
%
%       571.06   873.91   587.67   406.10   395.47   335.39   426.78   545.73   534.31   487.81   599.39   555.64   595.26   535.12   566.72   419.41
%
%     Columns 17 through 30:
%
%       413.43   581.49   608.56   544.89   354.99   343.75   478.57   709.40   543.41   623.79   450.16   326.30   532.02   514.18
%
%    weight =
%
%     Columns 1 through 29:
%
%       50   50   49   50   50   50   50   50   49   50   49   49   50   50   49   49   50   48   50   50   50   50   50   50   48   46   49   50   50
%
%     Column 30:
%
%       49
%
%    items =
%    {
%      [1,1] =
%
%          8.0000   12.4850
%         10.0000   40.3960
%          1.0000   20.3510
%          5.0000   70.4110
%          4.0000   69.8660
%          3.0000    7.9810
%          7.0000   44.0040
%          1.0000    8.3670
%          1.0000   92.7650
%          3.0000   10.5410
%          1.0000   95.0040
%          5.0000   21.9250
%          1.0000   76.9680
%
%      [1,2] =
%
%          1.0000   23.1200
%          3.0000   21.3490
%          1.0000   64.1410
%          4.0000   67.2780
%          1.0000   72.2280
%          5.0000   87.4870
%          2.0000   88.6910
%          5.0000   22.2780
%         10.0000   67.6650
%          4.0000   94.8750
%          3.0000   61.6310
%          7.0000   92.1200
%          3.0000   27.3810
%          1.0000   83.6680
%    
%
% ... and so on (rest of the output is omitted). The function loadInstances can be used to read instances from files (see help for more details)

function [solutionTotal] = KnapsackMainSet (itemsSet, maxWeightSet, rules, actions, featureIDs, varargin)  
  
  isTransformed = false; M = -1; W = -1;
  if ~isempty(varargin), isTransformed = true; M = varargin{1}; W = varargin{2}; end
  
  nbInstances = length(itemsSet);
  solutionTotal(nbInstances) = createEmptyKP(1);
%   for idx = 1 : nbInstances
%   parfor (idx = 1 : nbInstances, 0) % For debugging
  parfor idx = 1 : nbInstances
    kp = createEmptyKP(maxWeightSet(idx));
    if isTransformed        
      [solution]      = KnapsackMain(kp, itemsSet{idx}, rules, actions, featureIDs, M, W);
    else
      [solution]      = KnapsackMain(kp, itemsSet{idx}, rules, actions, featureIDs);
    end
    solutionTotal(idx) = solution;
  end
end