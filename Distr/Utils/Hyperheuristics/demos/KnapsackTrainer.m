% KnapsackTrainer
%% Basics
clc 
clear
close all

if ismac
    % Code to run on Mac platform
    addpath(genpath('../../../Utils'));
elseif isunix
    % Code to run on Linux platform
    addpath(genpath('../../../Utils'));
elseif ispc
    % Code to run on Windows platform
    addpath(genpath('..\..\..\Utils'));
else
    error('Platform not supported')
end

rng(0) % resets seed for debugging

%% Load instances
% // Load instances from file
% [wholeData, wholeMaxWeight] = loadInstances('../domains/knapsack/instances/kpInstances/Hard-Pisinger-200/');
% [wholeData, wholeMaxWeight] = loadInstances('../domains/knapsack/instances/kpInstances/ReducedSetForTesting/');
[wholeData, wholeMaxWeight] = loadInstances('../domains/knapsack/instances/kpInstances/AdvancedReducedEasy/');

% // Sort randomly the indices
scrambledIndex = randperm(length(wholeData));

% // Specify how many instances are used for training
% nbTrainInstances = 5;
% nbTrainInstances = 30;
nbTrainInstances = 80;

% // Switch for data using randperm
trainIndex = scrambledIndex(1:nbTrainInstances); 
testIndex = scrambledIndex(nbTrainInstances+1:end);

%% Set up runs

% // Define training parameters
nbRules = 2; nbActions = 4; 
% nbRules = 5; nbActions = 4;

featureIDs = [1 4 7];
% featureIDs = 1:7; % Original features
% featureIDs = 8:14; M = [0.25 0.25 0.1 0.6 0.15 0.1 0.3 0.35]; W = [0.2 0.2 0.1 0.2 0.15 0.1 0.15 0.2]; % Linear features
% featureIDs = 15:21; M = [0.25 0.25 0.1 0.6 0.15 0.1 0.3 0.35]; W = [0.2 0.2 0.1 0.2 0.15 0.1 0.15 0.2]; % S-shaped features
% featureIDs = 15:21; M = [0.25885   0.32000   0.28920   0.55940   0.16835   0.24490   0.29730   0.33845]; ...
%                     W = [0.16565   0.18860   0.21920   0.19560   0.16835   0.24490   0.20270   0.17175]; % S-shaped features with CEC(CIM) values
% featureIDs = 8:14; M = [0.25885   0.32000   0.28920   0.55940   0.16835   0.24490   0.29730   0.33845]; ...
%                     W = [0.16565   0.18860   0.21920   0.19560   0.16835   0.24490   0.20270   0.17175]; % Linear features with CEC(CIM) values

nbFeatures = length(featureIDs);

isTransformed = false; if any(featureIDs > 7), isTransformed = true; end

if isTransformed
  fObj = @(X)-KnapsackSelector(X, {wholeData{trainIndex}}, wholeMaxWeight(trainIndex), nbRules, nbActions, featureIDs, M, W); 
else
  fObj = @(X)-KnapsackSelector(X, {wholeData{trainIndex}}, wholeMaxWeight(trainIndex), nbRules, nbActions, featureIDs); 
end

%% Training Selector
% Training parameters

% - General parameters
parameters.verboseMode      = true; % Print mode of UPSO
parameters.maxIter          = 60; % Default: min(100,10*NUMBER_OF_VARIABLES).
% parameters.maxStagIter      = min(30,0.1*parameters.maxIter); % Default
parameters.maxStagIter      = parameters.maxIter/20; 
% parameters.randomSeed       = -1; % Default
% parameters.knownMin         = -inf; % Default
parameters.stagTol          = 1e-3;

% - Specific UPSO parameters
parameters.populationSize   = 20;
parameters.unifyFactor      = 0.2;
parameters.globalConf       = 2.833;
parameters.selfConf         = 2.833;
% parameters.phiFactor        = parameters.globalConf + parameters.selfConf;
% parameters.kappaFactor      = 1; 
% parameters.chiFactor        = 2*parameters.kappaFactor/...
%     abs(2-parameters.phiFactor-sqrt(parameters.phiFactor^2-4*parameters.phiFactor))

% - Selector parameters
parameters.numFeatures      = nbFeatures;
parameters.numRules         = nbRules; % +1 because of actions
parameters.Nrep             = 1; % Number of repetitions for training
parameters.printMode        = false;

% Trains selector
fprintf('Training over %d instances with %d rules...\n', length(trainIndex), nbRules);
[fP,P,details] = KnapsackTrainSelector(fObj, parameters); 
fprintf('Training complete in %d seconds!!!\n', details.elapsedTime);

% Prints best results of training
fP,P, 
% Prints results as set of rules (selector)
[R,A] = getRules(P,nbRules,nbActions); 
printRules(R,A)

%% Solves remaining instances with results of training
fprintf('Solving over %d instances (test set)...\n', length(testIndex));
totalFitness = KnapsackSelector(P, {wholeData{testIndex}}, wholeMaxWeight(testIndex), nbRules, nbActions, featureIDs);
fprintf('\n\n\t\t%.5e\n', totalFitness);

% % Solves all instances with results of training
% fprintf('Solving over all instances...\n');
% %totalFitness = KnapsackSelector(P, wholeData, wholeMaxWeight, nbRules, nbActions);
% if isTransformed
%   totalFitness = KnapsackSelector(P, wholeData2, wholeMaxWeight2, nbRules, nbActions, featureIDs, M, W); 
% else
%   totalFitness = KnapsackSelector(P, wholeData2, wholeMaxWeight2, nbRules, nbActions, featureIDs);
% end
% fprintf('\n\n\t\t%.5e\n', totalFitness);


%KnapsackSelector(P, {wholeData{1}}, wholeMaxWeight(1), nbRules, nbActions);