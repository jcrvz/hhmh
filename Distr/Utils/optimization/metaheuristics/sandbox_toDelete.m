% =====================================================================================
% INPUT ARGUMENTS:
% =====================================================================================
%
% * objectiveFunction   - The minimisation problem to be solved. It must be a script 
%                         or handle function. This parameter is mandatory.
%
% * popInitRange        - Simple boundaries of the search space. It must be a 
%                         NUMBER_OF_VARIABLES-times-2 matrix where the first and 
%                         second columns are the lower and upper boundaries per 
%                         dimension (row), respectively. If it is not provided, it 
%                         is defined as [-1,1;-1,1] assuming two dimensions, and the 
%                         par.simpleConstraints = false, for an unconstrained search.
%
% * par (optional):     - Parameter for running the algorithm.
%
%   General parameters
%   ------------------
%   verboseMode         - Additional details about what the algorithm is doing 
%                         are printed. It must be boolean. The default is true.
%   populationSize      - Number of agents conforming the population. It must be 
%                         an integer greater than 1. The default is 
%                         min(100,10*NUMBER_OF_VARIABLES).
%   simpleConstraints   - Define if the search space is constrained by the provided
%                         boundaries in popInitRange. It must be boolean. The default 
%                         is true.
%   randomSeed          - It is the seed for the random number generator. It must be 
%                         a real a non-negative integer. If a negative number is 
%                         provided, the default value, 'SHUFFLE', is chosen.
%   previousPopulation  - Allows an external population to be used as the initial 
%                         population. If previousPopulation is not provided, then the 
%                         method initialise the population. It must be a 
%                         NUMBER_OF_VARIABLES-times-POPULATION_SIZE matrix.
%   knownMin            - Provides the known minimum value for stopping the algorithm
%                         The default value is -INF.
%   maxIter             - Maximum number of iteratios. It must be a non-negative 
%                         integer. The default value is 100.
%   maxStagIter         - Maximum number of stagnated iterations. It must be a 
%                         non-negative integer. The default value is 10% of maxIter.
%   funcTol             - Function variation tolerance for the stopping criteria. 
%                         It must be positive. The default value is 1e-12.
%   solTol              - Solution variation tolerance for the stopping criteria. 
%                         It must be positive. The default value is 1e-16.
%   fitnessLimit        - Fitness threshold value (using the provided known minimum 
%                         value, if it is finite and real). It must be positive. 
%                         The default value is 1e-6.
%
%   PSO specific parameters
%   -----------------------
%   SelfConf            - PSO's self-confidence coefficient. It must be a real greater 
%                         than 0. The default value is 1.49.
%   GlobalConf          - PSO's global-confidence coefficient. It must be a real greater 
%                         than 0. The default value is 1.49.
%   kappaFactor         - PSO's factor for constriction factor calculation. It must be 
%                         a real greater than 0. The default value is 1.
%   phiFactor           - PSO's factor for constriction factor calculation. It must be 
%                         a real greater than 0. The default value is given by
%                         SELF_CONFIDENCE + GLOBAL_CONFIDENCE.
%   chiFactor           - PSO's constriction factor calculation. It must be a real 
%                         greater than 0. The default value is given by
%                         2*KAPPA / ABS ( 2 - PHI - SQRT( PHI^2 - 4*PHI ) ).
%   initialVelMode      - Defines how velocities are initialised. It must be a string
%                         such as 'rand' or 'zero'. The default is 'rand'.
%
%   
%   UPSO specific parameters
%   ------------------------
%   unifyFactor         - UPSO's Unification factor. It must be a real between [0,1].
%                         The default value is 0.5.
%   neighbourhood       - Population neighbourhood topology. It must be a binary
%                         POPULATION_SIZE-times-POPULATION_SIZE matrix. The default 
%                         value is the matrix of a (3-neighbours) ring topology.
%
% =====================================================================================
% OUTPUT ARGUMENTS:
% =====================================================================================
%
% * globalBestPosition  - Global best position is a NUMBER_OF_VARIABLES-vector 
%                         containing the latest best solution found.
% * globalBestFitness   - Global best fitness is the objectiveFunction evaluated in 
%                         globalBestPosition.
% * details:            - Details of the performed procedure.
%
%   Fields of details
%   -----------------
%   elapsedTime         - Overall procedure elapsed time in seconds.
%   functionEvaluations - Total performed evaluations number
%   performedIterations - Total iteration number
%   procedureEvolution: - A structure with information about the algorithm
%                         performance per interation.
%
%       Fields of procedureEvolution 
%       ----------------------------
%       solution        - A NUMBER_OF_VARIABLES-times-NUMBER_OF_ITERATIONS matrix 
%                         which contains the best position found per iteration.
%       fitness         - Structure with information about the fitness improvement per 
%                         iteration.
%
%           Fields of fitness 
%           -----------------
%           raw         - A NUMBER_OF_ITERATIONS vector which contains the best fitness 
%                         values found per iteration.
%           mean        - A NUMBER_OF_ITERATIONS vector which contains the average value
%                         of best fitness values found until the corresponding iteration.
%           stdv        - A NUMBER_OF_ITERATIONS vector which contains the standard 
%                         deviation value of best fitness values found until the 
%                         corresponding iteration.
%
%   stoppingFlags:      - A structure with the latest values of the stopping flags of 
%                         criteria
%
%       Fields of stoppingFlags 
%       -----------------------
%       iterationFlag   - Criterion 1: Has iteration counter reached the maximum number?
%       maxStagIter     - Criterion 2: Is the procedure stagnated?
%       fitnessLimitFlag- Criterion 3: Has the best fitness reached the fitness limit?
%

Na = 5;
SSw     = diag(ones(Na,1)) + diag(ones(Na - 1,1),1) + ...
    diag(ones(Na - 1,1),-1) + diag(ones(1),Na - 1) + ...
    diag(ones(1),(1 - Na));
SSw(SSw==0)=nan

fPi = [1 20 -3 100 3.3]'


repmat(fPi,1,Na).*SSw

[mV,mP] = min(repmat(fPi,1,Na).*SSw)