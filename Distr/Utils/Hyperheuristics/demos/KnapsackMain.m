% Function for the Knapsack (KP) problem (solves a single instance)
% -------------------------------------------------
% Inputs
% -------------------------------------------------
% knapsack: Empty knapsack created with createEmptyKP(...)
% items:    Data matrix. Row: Item. Columns: Weight, Profit
% rules:    Matrix containing the set of rules. One row per rule and one column per feature
% actions:  Column vector with action indices, starting at one
% featureIDs: IDs of features to use
%
% -------------------------------------------------
% Outputs
% -------------------------------------------------
% solution: Structure containing the solution. It has the following fields: fitness, weight, items.
%           The first two give the total fitness and total weight. The last one gives the packed items.
%
% -------------------------------------------------
% Example
% -------------------------------------------------
% ti = [3 10;4 2;15 1;5 20];
% mw = 20;
% ts = [0.5*ones(1,7)];
% ta = [2];
% KnapsackMain(ti,mw,ts,ta)
%
% Returns:
% scalar structure containing the fields:
%
%    fitness =  32
%    weight =  12
%    items =
%
%            5   20
%            3   10
%            4    2
%
% The function loadInstances can be used to read instances from files (see help for more details)

function [knapsack] = KnapsackMain(knapsack, items, rules, actions, featureIDs, varargin)
% addpath(genpath('../../../Utils'));
isTransformed = false; M = -1; W = -1;
if ~isempty(varargin), isTransformed = true; M = varargin{1}; W = varargin{2}; end
items_ = items;
% [nbRules, nbFeatures] = size(rules);
while true
    % Extract features
    if isTransformed
%         curr_features = getKPFeatures(items_, featureIDs, M, W); % Asks for the desired features
        curr_features = getKPFeatures(items_, featureIDs, M, W, knapsack); % Temp... just to test if solution yields good data
    else
        curr_features = getKPFeatures(items_, featureIDs); % Asks for the desired features
    end
    % Compare to rules
    rule_dist = measureRules(curr_features, rules);
    % Decide action
    [~, action_id] = min(rule_dist);
    action = actions(action_id);
    % Apply action
    itemID = selectItemKP(items_, action);
    [knapsack, items_, exFlag] = storeKPItem(knapsack, items_, itemID);
%     exFlag;
    % Check if item did not fit
    if exFlag == -1
        items_ = [items_(1:itemID-1,:);...
                    items_(itemID+1:end,:)];
    end
    % Check convergence
    if isempty(items_)
        break;
    end
    %fflush(stdout);
end
knapsack = knapsackFitness(knapsack);