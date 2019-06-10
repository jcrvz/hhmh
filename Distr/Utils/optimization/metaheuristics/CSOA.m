function [bestAgent,bestFitness,details] = CSOA(ObjectiveFunction,PopInitRange)
%
%% Initialise DiscoveryProbrameters
% - Simple constraints activation
applySimpleConst    = true;

% - Number of Variables
numberOfVariables   = size(PopInitRange,1);

% - Number of agents or DiscoveryProbrticles
PopulationSize      = min(100,10*numberOfVariables);

% - Cuckoo Search's DiscoveryProbrameters
Beta = 1.5; % Beta DiscoveryProbrameter to calculate Sigma
DiscoveryProb = 0.75; % Probability of discover eggs
ScaleFactor = 0.1; % Step size scale factor

% - Stopping criteria limits
MaxIterations       = 500*numberOfVariables;
MaxStagIterations   = 50;
% StepTolerance       = 1e-6;
FunctionTolerance   = 1e-6;
FitnessLimit        = eps;

%% Initialise the algorithm's variables and DiscoveryProbrameters
% Distribute the initial population uniformly along Search Space
[CurrentPopulation,PopulationSize] = InitialisePopulation('rand');

% Adjust the boundaries of the optimisation's domain in matrix form
LowerBound = repmat(PopInitRange(:,1),1,PopulationSize);
UpperBound = repmat(PopInitRange(:,2),1,PopulationSize);

% - Calculate std dev of u (Mantegna's algorithm)
Sigma = (gamma(1 + Beta)*sin(pi*Beta/2)/(gamma((1 + Beta)/2)*Beta*...
    2^((Beta - 1)/2)))^(1/Beta);

% Calculate the initial fitness values
CurrentEvaluatedFunction = EvalFunction();

% Found the best position to be rotation centre
[bestFitness,bestFitnessId] = min(CurrentEvaluatedFunction);
bestAgent = CurrentPopulation(:,bestFitnessId);

% Store the previous best solution
previousBestFitness = inf;
% previousBestAgent   = inf(size(bestAgent));

% Store the previous population
previousPopulation = CurrentPopulation;
previousEvaluatedFunction = CurrentEvaluatedFunction;

% Save the current best fitness
historicalFitness    = nan(1,MaxIterations+1);
historicalFitness(1) = bestFitness;

% Set the iteration and stagnation counters
Iterations = 0;
% StagIterations = 0;
bestFvalsWindow = nan(MaxStagIterations, 1);

% Set the criteria flag as false
criteria = false;

%% Main proceduce: Cuckoo Search Optimisation Algorithm (CSOA)
% - Start a time counter
timerVal = tic;

% - Repeat the following process till criteria becomes true
while ~criteria
    
    % Update the population's locations
    LevyFlightUpdating();
    
    % Check if the DiscoveryProbrticle is inside the search
    SimpleConstraints(applySimpleConst);
    
    % Evaluate objective function in new positions
    CurrentEvaluatedFunction = EvalFunction();
    
    % Update population to preserve improvements
    %UpdatePopulation();
    
    % Update the population's locations
    DiscoveryUpdating();
    
    % Check if the DiscoveryProbrticle is inside the search 
    SimpleConstraints(applySimpleConst);
    
    % Evaluate objective function in new positions
    CurrentEvaluatedFunction = EvalFunction();
    
    % Update population to preserve improvements
    UpdatePopulation();
        
    % Found the best position to be rotation center
    [bestFitness,bestFitnessId] = min(CurrentEvaluatedFunction);
    
    % Is it better than the previous solution?
    if bestFitness < previousBestFitness
        bestAgent = CurrentPopulation(:,bestFitnessId);     
    else
        bestFitness = previousBestFitness;
    end 
    
    bestFvalsWindow(1+mod(Iterations-1,MaxStagIterations)) = bestFitness;
    
    % Update the criteria flag
    [criteria,stoppingFlags] = UpdateCriteria();
    
    % Store the previous best solution
    previousBestFitness = bestFitness;
%     previousBestAgent   = bestAgent;
    
    % Save the current best fitness
    historicalFitness(Iterations) = bestFitness;
    
    % <TO-DELETE> Plot
%     Color = exp(-Iterations/50) - 1e-3;
%     plot(CurrentPopulation(1,:),CurrentPopulation(2,:),'.', ...
%        'Color',Color*[0.5 0.5 1]); hold on;
%     plot(bestAgent(1),bestAgent(2),'*','Color',Color*[1 0 0]); %hold off;
%     axis([0 1 0 1]); pause(0.1); 
%     getframe(gcf);
    % </TO-DELETE>
    
end

% - End the time counter
elapsedTime = toc(timerVal);

% - Rescale the bestAgent to the problem's domain
bestAgent = RescaleVariables(bestAgent);

% - Store details about the performed procedure
details = struct('elapsedTime',elapsedTime, ...
    'functionEvaluations',Iterations*(PopulationSize + 1), ...
    'performedIterations',Iterations, ...
    'historicalFitness',historicalFitness(~isnan(historicalFitness)));

for fn = fieldnames(stoppingFlags)'
   details.(fn{1}) = stoppingFlags.(fn{1});
end

%% Definition of functions used by this algorithm

% ------------------------------------------------------------------------
% Function for the initial population distribution
% ------------------------------------------------------------------------

    function [InitialPopulation,NewPop] = InitialisePopulation(kindOf)
        % Calculate the initial positions, y (in [0,1])
        if nargin < 1, kindOf = 'rand'; end
        switch kindOf
            case 'grid'
                % Find points per dimensions and correct the population (if required)
                [InitialPopulation,NewPop] = ...
                    gridDistribution(numberOfVariables,PopulationSize);
                
            otherwise
                NewPop = PopulationSize;
                InitialPopulation = rand(numberOfVariables,PopulationSize);
        end
    end

% ------------------------------------------------------------------------
% Function for rescaling the variables
% ------------------------------------------------------------------------

    function [RescaledVar] = RescaleVariables(Vars)
        % Create the Synthesis equation to transform y (in [0,1]) to x 
        %   (in [lower,upper]). Both x and y are numOfVar--times--PopSize
        
        % Define the boundaries for each dimension
        %PopInitRange    = [min(PopInitRange,[],2) max(PopInitRange,[],2)];
        
        % Variables can be one agent or the entire population, then the
        % function returns one vector or a matrix.
        if size(Vars,2) == PopulationSize                     
            RescaledVar = LowerBound + Vars.*(UpperBound - LowerBound);
        else
            RescaledVar = LowerBound(:,1) + ...
                Vars.*(UpperBound(:,1) - LowerBound(:,1));
        end
    end

% ------------------------------------------------------------------------
% Function for the objective function evaluation
% ------------------------------------------------------------------------

    function EvaluatedFunction = EvalFunction()
        % Check if it is the first evaluation
        if exist('EvaluatedFunction','var') == false
            EvaluatedFunction = nan(PopulationSize,1);
        end
        
        % Rescale the population values to the problem's domain
        rescaledPopulation = RescaleVariables(CurrentPopulation);
        
        % Evaluate each agent's position into the objective function
        for AgentId = 1 : PopulationSize
            EvaluatedFunction(AgentId) = ...
                ObjectiveFunction(rescaledPopulation(:,AgentId));
        end
    end

% ------------------------------------------------------------------------
% Function for the Levy's flight updating
% ------------------------------------------------------------------------

    function LevyFlightUpdating()
        % Generate random numbers with the uniform distribution
        valsU = randn(numberOfVariables,PopulationSize);
        
        % Generate random numbers with the Normal distribution
        valsN1 = rand(numberOfVariables,PopulationSize)*Sigma;
        valsN2 = abs(randn(numberOfVariables,PopulationSize));
        
        % Obtain random numbers with the Levy distribution
        LevyRand = valsN1./(valsN2.^(1/Beta));
        
        % Calculate the Levy's flight for the entire population
        CurrentPopulation = CurrentPopulation + ...
            ScaleFactor*valsU.*LevyRand.*...
            (CurrentPopulation - repmat(bestAgent,1,PopulationSize));        
    end

% ------------------------------------------------------------------------
% Function for the Discovery updating
% ------------------------------------------------------------------------

    function DiscoveryUpdating()
        % Generate random numbers with the uniform distribution
        valsU1 = rand(numberOfVariables,PopulationSize);
        valsU2 = rand(numberOfVariables,PopulationSize);
        
        % Update population checking if each pal was discovered
        CurrentPopulation = CurrentPopulation + valsU1.*...
            (CurrentPopulation(:,randperm(PopulationSize)) - ...
            CurrentPopulation(:,randperm(PopulationSize))).* ...
            (valsU2 > DiscoveryProb);
    end

% ------------------------------------------------------------------------
% Function for the simple constraints verification
% ------------------------------------------------------------------------

    function SimpleConstraints(activation)
        if activation == true
            CurrentPopulation(CurrentPopulation < 0) = 0;
            CurrentPopulation(CurrentPopulation > 1) = 1;
        end
    end

% ------------------------------------------------------------------------
% Function for the update population
% ------------------------------------------------------------------------

    function UpdatePopulation()
        % Is the population improving?
        improvingCond = ones(numberOfVariables,1)*...
            (CurrentEvaluatedFunction < previousEvaluatedFunction)';
        
        % Update the population according with its goodness 
        CurrentPopulation(~improvingCond) = previousPopulation(~improvingCond);
        CurrentEvaluatedFunction = min(CurrentEvaluatedFunction, ...
            previousEvaluatedFunction);
        
        % Save the best values
        previousPopulation = CurrentPopulation;
        previousEvaluatedFunction = CurrentEvaluatedFunction;
    end

% ------------------------------------------------------------------------
% Function for the criteria evaluation
% ------------------------------------------------------------------------

    function [criteria,stoppingFlags] = UpdateCriteria()
        % Criterion 1: Is the iteration counter reached the maximum number?
        Iterations     = Iterations + 1;
        stoppingFlags.IterFlag = Iterations >= MaxIterations;
        
        % Criterion 2: Is the procedure stagnated?
        % -> 2.1 Equal fitness
%         equalFitness = (globalBestFitness == previousBestFitness);
%         % -> 2.2 Tiny difference between current and previous solution
%         nonZeroStepTol = 1 + (norm(previousBestAgent) < 1e-6)*norm(previousBestAgent);
%         closeSolution = (norm(previousBestAgent - globalBestPosition) < ...
%             StepTolerance*nonZeroStepTol);
%         % -> 2.3 Tiny difference between current and previous fitness
%         nonZeroFuncTol = 1 + (abs(previousBestFitness) < 1e-3)*abs(previousBestFitness);
%         closeFitness = (abs(previousBestFitness - globalBestFitness) < ...
%             FunctionTolerance*nonZeroFuncTol);
%         if equalFitness || closeSolution || closeFitness
%             StagIterations = StagIterations + 1;
%         end
%         stoppingFlags.StagFlag = (StagIterations >= MaxStagIterations);
        
        iterationIndex = 1+mod(Iterations-1,MaxStagIterations);
        if Iterations > MaxStagIterations
            if iterationIndex == MaxStagIterations
                % The window runs from index 1:iterationIndex
                maxBestFvalsWindow = bestFvalsWindow(1);
            else
                % The window runs from [iterationIndex+1:end, 1:iterationIndex]
                maxBestFvalsWindow = bestFvalsWindow(iterationIndex+1);
            end
            funChange = abs(maxBestFvalsWindow-bestFitness)/max(1,abs(bestFitness));
        else
            funChange = Inf;
        end
        stoppingFlags.StagFlag = (funChange < FunctionTolerance);
        
        %         bestFvalsWindow(1+mod(state.Iteration-1,options.StallIterLimit)) = min(state.IndividualBestFvals);
        
        % Criterion 5: Has the best fitness reached the fitness limit?
        stoppingFlags.FitnessLimitFlag = abs(bestFitness) < FitnessLimit;
        
        % Summary of criteria
        criteria = stoppingFlags.IterFlag || stoppingFlags.StagFlag || ...
            stoppingFlags.FitnessLimitFlag;
    end

end
