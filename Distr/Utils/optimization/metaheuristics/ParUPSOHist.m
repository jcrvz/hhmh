%% MODIFIED UPSO ALGORITHM - ALGORITMO MODIFICADO DE UPSO
%
% Developed by Jorge Cruz, Eng.
% Advised by Rodrigo Correa, Ph.D. and Iván Amaya, Ph.D(c).
%
% CEMOS Research Group - Grupo de Investigación CEMOS
% Universidad Industrial de Santander%
% Bucaramanga - Santander - Colombia
% 14 - Oct - 2012
%
% 2014-apr ver :: 1

function [Pg,fPg,details, varargout] = ParUPSOHist(fObj,bnd,parameters,oldX)
% Check which call was made
hasInner = false;
if nargout > 3, hasInner = true; end 

% Read dimensions
Nd      = size(bnd,1);

% Initial positions 
bnd     = [min(bnd,[],2) max(bnd,[],2)];
bli(1,:) = bnd(:,1); bdif(1,:) = bnd(:,2)' - bli;

% Read parameters
if nargin < 3
    Na      = 100;
    cp      = 2.0;
    cg      = 2.5;
    chi     = 0.6;
    u       = 0.5;
    
    Tol     = 1;
    M       = 10000;
    msat    = 100;
    
    seed    = -1;
    rand("state");
    
    knownMin = 0;
    
    X       = repmat(bli,Na,1) + rand(Na,Nd).*repmat(bdif,Na,1);
    
    printMode = true;
	toHistory = false;
else
    Na      = parameters.NAG;
    cp      = parameters.CP;
    cg      = parameters.CG;
    chi     = parameters.CHI;
    u       = parameters.U;
    
    Tol     = parameters.TOL;
    M       = parameters.MITE;
    msat    = parameters.MSAT;
    
    knownMin = parameters.knownMin;
    
    seed    = parameters.SEED;
    if seed == -1, rng('shuffle'); else, rng(seed); end
    
    if nargin == 4    
      X = oldX;
    else
      X       = repmat(bli,Na,1) + rand(Na,Nd).*repmat(bdif,Na,1);
    end
    
    printMode = parameters.printMode;
	toHistory = parameters.toHistory;
end

%% -- Debug
%fprintf("\nStarting Pop:");
%      for idp = 1 : parameters.NAG
%        fprintf("\t%.4f",X(idp,:));
%        fprintf("\n");
%      end
%fflush(stdout);

if toHistory, histGb = zeros(M+1,Nd); end


% Pre-allocate some variables
fX      = zeros(Na,1);
U       = zeros(Na,Nd);
if hasInner, allInternals(Na) = struct('solution',-1,'details',-1); end

% Define the ring topology
SSw     = diag(ones(Na,1)) + diag(ones(Na - 1,1),1) + ...
    diag(ones(Na - 1,1),-1) + diag(ones(1),Na - 1) + ...
    diag(ones(1),(1 - Na));


% Initial values
% for i = 1 : Na 
parfor i = 1 : Na
	if hasInner
		[fX(i,1), internalSol, internalDetails] = fObj(X(i,:)); 
		allInternals(i) = struct('solution',internalSol,'details',internalDetails);
	else
		fX(i,1) = fObj(X(i,:)); 
	end
	
end

fPi     = fX; Pi = X;

% if Na == 5, 
%   pause(0.1)
% end


% Found the neighbourhood best position (initial step)
[~,gi]  = min(repmat(fPi,1,Na).*(SSw == 1) + ~(SSw == 1)*100);
Pgi     = Pi(gi,:);

% Found the swarm best position (initial step)
[fPg,g] = min(fPi); Pg = X(g,:);

if toHistory, histGb(1,:) = Pg; end

% Save initial solution
if hasInner
    varargout{1} = allInternals(g); 
end


% % For debugging parfor
% details = struct('time',[],'fevs',[],'steps',[],...
%     'outmsg',[],'favg',[],'fstd',[],...
%     'conc',[], 'spread',[], 'disper',[], 'pop',[],...
%     'a', []);
% return


% Set auxiliar variables
steps    = 0;
msatc   = 0;
sumAVG  = 0;
sumSD   = 0;

% Plot objective function value evolution
% topl = [];

% Print best fitness after initialization
if printMode, fprintf('\tBest fitness at startup: %.4e\n', fPg), end
% fflush(stdout);




% ------------------------------------
concentration = GetConcentration(X, Pg, (bnd(:,2)-bnd(:,1))');
conc   = TransformConcentration(concentration, 0);
spread = TransformConcentration(concentration, 1);
disper = TransformConcentration(concentration, 2);
%fprintf("%.5f\t",concentration, conc, spread, disper);
%fprintf("\n");
%fflush(stdout);
% ------------------------------------


bndL = bnd(:,1); 
bndU = bnd(:,2);

%% Main process
tic,
while steps < M && msatc < msat && fPg > knownMin
    % Update step
    steps    = steps + 1;
    if printMode, fprintf('\tSteps = %d ...\tFitness = %.5e ...\n', steps, fPg); end
%     fflush(stdout);
    
    % Update the global and local velocities for each particle
    Gu      = chi*(U + cp*rand(Na,Nd).*(Pi - X) + ...
        cg*rand(Na,Nd).*(repmat(Pg,Na,1) - X));
    Lu      = chi*(U + cp*rand(Na,Nd).*(Pi - X) + ...
        cg*rand(Na,Nd).*(Pgi - X));
    
    % Update the total velocity for each particle
    U       = (1 - u)*Lu + Gu*u;
    
    % Update the position for each particle
    X       = X + U;
    
%     for i = 1 : Na
    parfor i = 1 : Na
%     parfor (i = 1 : Na,0)
        
        % Check if the particle is in search space
        newX = X(i,:);
        
        % Another way
%         targetIndex = newX < bndL' | newX > bndU';
%         newX(targetIndex)   = bli(targetIndex) + rand(1,sum(targetIndex)).*bdif(targetIndex);
       
        for d = 1 : Nd
            if newX(d) < bndL(d) || newX(d) > bndU(d)
                newX(d)   = bli(d) + rand*bdif(d);
            end
        end
        X(i,:) = newX;
        
        % Evaluate objective function in the new position
        if hasInner
            [fX(i), internalSol, internalDetails] = fObj(X(i,:));
            allInternals(i).solution = internalSol;
            allInternals(i).details = internalDetails;
        else
            fX(i)   = fObj(X(i,:));
        end
    end
    
    % Found the best position for each particle
    cond    = fX < fPi;
    Pi      = repmat(cond,1,Nd).*X + repmat(~cond,1,Nd).*Pi;
    fPi     = min(fX,fPi);
    
    % Found the neighbourhood's best position
    [~,gi]  = min(repmat(fPi,1,Na).*(SSw == 1) + ~(SSw == 1)*100);
    Pgi     = Pi(gi,:);
    
    fPg_ = fPg;
    % Found the swarm best position
    [fPg,g] = min(fPi);
    if fPg < fPg_
        Pg = X(g,:);
        if hasInner
            varargout{1} = allInternals(g);
        end
    else
        fPg = fPg_;
    end
	
	% Updates historical values 
	if toHistory, histGb(steps+1,:) = Pg; end
    
    % Statistical
    sumAVG      = sumAVG + fPg;
    sumSD       = sumSD + fPg^2;
    currAVG     = sumAVG/steps;
    currSD      = sqrt(sumSD/steps - currAVG^2);    
    
    % Stop criteria
    if abs(fPg) > abs(currAVG) - Tol*currSD && fPg == fPg_,
        msatc = msatc + 1;
    else
        msatc = 0;
    end
    
%     topl        = [topl; [steps,fPg,currAVG,currSD]];
end
t       = toc;

% Plot objective function value evolution
% figure('Name','Enjambre','Color','White');
% hold on, errorbar(topl(:,1),topl(:,3),topl(:,4),'go',...
%     'MarkerFaceColor','b','MarkerEdgeColor','b');
% plot(topl(:,1),topl(:,2),'k-.','LineWidth',2); hold off,
% xlim([1 steps]), legend('fAVG','fBest'), %getframe(gcf);

if steps < M,   outmsg = 1;
else            outmsg = 0;
end



% ------------------------------------
concentration = GetConcentration(X, Pg, (bnd(:,2)-bnd(:,1))');
conc   = TransformConcentration(concentration, 0);
spread = TransformConcentration(concentration, 1);
disper = TransformConcentration(concentration, 2);
% ------------------------------------



details = struct('time',t,'fevs',(steps+1)*Na,'steps',steps,...
    'outmsg',outmsg,'favg',currAVG,'fstd',currSD,...
    'conc',conc, 'spread',spread, 'disper',disper, 'pop',Pi,...
    'a', concentration);

if toHistory, details.bestHistory = histGb(1:steps,:); end 
