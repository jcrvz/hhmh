function [globalBestPosition,globalBestFitness,details] = ...
    PSO(ObjectiveFunction, PopInitRange, Parameters)
if nargin < 3
    Parameters = struct();
    if nargin < 2
        Parameters.applySimpleConst = false;
        PopInitRange                = [-1,1;-1,1];
        if nargin < 1
            error('ObjectiveFunction must be provided!');
        end
    end
end

%% Initialise Parameters

% Number of Variables
Parameters.numberOfVariables = size(PopInitRange,1);

% Simple constraints activation
if ~isfield(Parameters,'applySimpleConst'), Parameters.applySimpleConst = true; end

% Number of agents
if ~isfield(Parameters,'PopulationSize'), ...
        Parameters.PopulationSize = min(100,10*Parameters.numberOfVariables); end

% Self-confidence coefficient
if ~isfield(Parameters,'SelfConfidence'), Parameters.SelfConfidence = 1.49; end

% Global-confidence coeffcient
if ~isfield(Parameters,'GlobalConfidence'), Parameters.GlobalConfidence = 1.49; end

% Parameters for the Constriction factor calculation
if ~isfield(Parameters,'KappaFactor'), Parameters.KappaFactor = 1; end

% Stopping criteria limits
if ~isfield(Parameters,'MaxIterations'), Parameters.MaxIterations = 100; end
if ~isfield(Parameters,'MaxStagIterations'), Parameters.MaxStagIterations = 101; end
if ~isfield(Parameters,'FunctionTolerance'), Parameters.FunctionTolerance = 1e-12; end
if ~isfield(Parameters,'FitnessLimit'), Parameters.FitnessLimit = 1e-16; end

%% Initialise the algorithm's variables and DiscoveryProbrameters

% Calculate the constriction factor
constPhi            = Parameters.SelfConfidence + Parameters.GlobalConfidence;
ConstrictionConstant= 2*Parameters.KappaFactor/abs(2 - constPhi - ...
    sqrt(constPhi^2 - 4*constPhi)); % Constriction constant (chi)

% Distribute the initial population uniformly along Search Space
[CurrentPopulation,CurrentVelocity,Parameters.PopulationSize] = InitialisePopulation('rand');

% Adjust the boundaries of the optimisation's domain in matrix form
LowerBound          = repmat(PopInitRange(:,1),1,Parameters.PopulationSize);
UpperBound          = repmat(PopInitRange(:,2),1,Parameters.PopulationSize);

% Calculate the initial fitness values
CurrentEvaluatedFunction    = EvalFunction();

% Store the particular best positions
particularBestPositions     = CurrentPopulation;
particularBestFitness       = CurrentEvaluatedFunction;

% Find the best position
[globalBestFitness,globalBestFitnessId] = min(particularBestFitness);
globalBestPosition  = particularBestPositions(:,globalBestFitnessId);

% Save the current best fitness
historicalFitness   = nan(1,Parameters.MaxIterations+1);
historicalFitness(1)= globalBestFitness;

% Set the iteration and stagnation counters
Iterations          = 0;
bestFvalsWindow     = nan(Parameters.MaxStagIterations, 1);

% Set the criteria flag as false
criteria = false;

%% Main proceduce: Unified Particle Swarm Optimisation (UPSO)
% - Start a time counter
timerVal = tic;

% - Repeat the following process till criteria becomes true
while ~criteria
    
    % Update the population's locations
    ParticleSwarmUpdating();
    
    % Check if the DiscoveryProbrticle is inside the search sDiscoveryProbce
    SimpleConstraints(Parameters.applySimpleConst);
    
    % Evaluate objective function in new positions
    CurrentEvaluatedFunction = EvalFunction();
    
    % Update the particular best positions
    improvingCond = logical(ones(Parameters.numberOfVariables,1)* ...
        (CurrentEvaluatedFunction < particularBestFitness)');
    particularBestPositions(improvingCond) = CurrentPopulation(improvingCond);
    particularBestFitness = min(CurrentEvaluatedFunction, ...
        particularBestFitness);
    
    % Store the previous best solution
    previousBestFitness = globalBestFitness;
    
    % Update the global best position
    [globalBestFitness,bestFitnessId] = min(particularBestFitness);
    if globalBestFitness < previousBestFitness
        globalBestPosition = particularBestPositions(:,bestFitnessId);
    else
        globalBestFitness = previousBestFitness;
    end
    
    bestFvalsWindow(1+mod(Iterations-1,Parameters.MaxStagIterations)) = globalBestFitness;
    
    % Update the criteria flag
    [criteria,stoppingFlags] = UpdateCriteria();
    
    
    % Save the current best fitness
    historicalFitness(Iterations) = globalBestFitness;
    
    % <TO-DELETE> Plot
    %     Color = exp(-Iterations/50) - 1e-2;
    %     plot(CurrentPopulation(1,:),CurrentPopulation(2,:),'.', ...
    %         'Color',Color*[0.5 0.5 1]); hold on;
    %     plot(globalBestPosition(1),globalBestPosition(2),'*','Color',Color*[1 0 0]); %hold off;
    %     axis([0 1 0 1]); pause(0.1);
    %     getframe(gcf);
    % </TO-DELETE>
    
end

% - End the time counter
elapsedTime = toc(timerVal);

% - Rescale the globalBestPosition to the problem's domain
globalBestPosition = RescaleVariables(globalBestPosition);

% - Store details about the performed procedure
details = struct('elapsedTime',elapsedTime, ...
    'functionEvaluations',Iterations*(Parameters.PopulationSize + 1), ...
    'performedIterations',Iterations, ...
    'historicalFitness',historicalFitness);%(~isnan(historicalFitness)));

for fn = fieldnames(stoppingFlags)'
    details.(fn{1}) = stoppingFlags.(fn{1});
end

%% Definition of functions used by this algorithm

% ------------------------------------------------------------------------
% Function for the initial population distribution
% ------------------------------------------------------------------------

    function [InitialPopulation,InitialVelocity,NewPop] = ...
            InitialisePopulation(kindOf)
        
        % Calculate the initial positions, y (in [0,1])
        if nargin < 1, kindOf = 'rand'; end
        switch kindOf
            case 'grid'
                % Find points per dimensions and fix the population (if so)
                [InitialPopulation,NewPop] = ...
                    gridDistribution(Parameters.numberOfVariables,Parameters.PopulationSize);
                
            otherwise
                NewPop = Parameters.PopulationSize;
                InitialPopulation = rand(Parameters.numberOfVariables,Parameters.PopulationSize);
        end
        
        MaxVelocity = (PopInitRange(:,2) - PopInitRange(:,1));
        InitialVelocity   = repmat(-MaxVelocity,1,Parameters.PopulationSize) + ...
            repmat(2*MaxVelocity,1,Parameters.PopulationSize) .* ...
            rand(Parameters.numberOfVariables,Parameters.PopulationSize);
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
        if size(Vars,2) == Parameters.PopulationSize
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
            EvaluatedFunction = nan(Parameters.PopulationSize,1);
        end
        
        % Rescale the population values to the problem's domain
        rescaledPopulation = RescaleVariables(CurrentPopulation);
        
        % Evaluate each agent's position into the objective function
        for AgentId = 1 : Parameters.PopulationSize
            EvaluatedFunction(AgentId) = ...
                ObjectiveFunction(rescaledPopulation(:,AgentId));
        end
    end

% ------------------------------------------------------------------------
% Function for the particle swarm updating
% ------------------------------------------------------------------------

    function ParticleSwarmUpdating()
        % Generate random numbers with the uniform distribution
        valsCp = Parameters.SelfConfidence*rand(Parameters.numberOfVariables,Parameters.PopulationSize);
        valsCg = Parameters.GlobalConfidence*rand(Parameters.numberOfVariables,Parameters.PopulationSize);        
        
        % Update the current velocity for each particle
        newVelocities = ConstrictionConstant*(CurrentVelocity + ...
            valsCp.*(particularBestPositions - CurrentPopulation) + ...
            valsCg.*(repmat(globalBestPosition,1,Parameters.PopulationSize) - ...
            CurrentPopulation));
        
        tfValid = all(isfinite(newVelocities), 1);
        CurrentVelocity(:,tfValid) = newVelocities(:,tfValid);
        
        % Update the position for each particle
        newPositions = CurrentPopulation + CurrentVelocity;
        
        tfValid = isfinite(newPositions);
        CurrentPopulation(tfValid) = newPositions(tfValid);
        
    end

% ------------------------------------------------------------------------
% Function for the simple constraints verification
% ------------------------------------------------------------------------

    function SimpleConstraints(activation)
        if activation == true
            CurrentPopulation(CurrentPopulation < 0) = 0;
            CurrentVelocity(CurrentPopulation < 0) = 0;
            
            CurrentPopulation(CurrentPopulation > 1) = 1;
            CurrentVelocity(CurrentPopulation > 1) = 0;
        end
    end

% ------------------------------------------------------------------------
% Function for the criteria evaluation
% ------------------------------------------------------------------------

    function [criteria,stoppingFlags] = UpdateCriteria()
        
        % Criterion 1: Is the iteration counter reached the maximum number?
        Iterations     = Iterations + 1;
        stoppingFlags.IterFlag = Iterations >= Parameters.MaxIterations;     
        
        % Criterion 2: Is the problem solution procedure stagnated?
%         iterationIndex = 1+mod(Iterations-1,Parameters.MaxStagIterations);
%         if Iterations > Parameters.MaxStagIterations
%             if iterationIndex == Parameters.MaxStagIterations
%                 % The window runs from index 1:iterationIndex
%                 maxBestFvalsWindow = bestFvalsWindow(1);
%             else
%                 % The window runs from [iterationIndex+1:end, 1:iterationIndex]
%                 maxBestFvalsWindow = bestFvalsWindow(iterationIndex+1);
%             end
%             funChange = abs(maxBestFvalsWindow-globalBestFitness)/max(1,abs(globalBestFitness));
%         else
%             funChange = Inf;
%         end
%         stoppingFlags.StagFlag = (funChange < Parameters.FunctionTolerance);
        
        % Criterion 3: Has the best fitness reached the fitness limit?
%         stoppingFlags.FitnessLimitFlag = abs(globalBestFitness) < Parameters.FitnessLimit;
        
        % Summary of criteria
        criteria = stoppingFlags.IterFlag;% || stoppingFlags.StagFlag || ...
%             stoppingFlags.FitnessLimitFlag;
    end

end
