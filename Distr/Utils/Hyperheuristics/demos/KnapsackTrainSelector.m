function [fPg, Pg, details] = KnapsackTrainSelector(fObj, parameters)
% --- Base parameters
if nargin < 2
    parameters.NAG   = 20;
    parameters.U     = 0.2;
    parameters.CG    = 2.833;
    parameters.CP    = 2.833;
    phi              = parameters.CP + parameters.CG;
    parameters.CHI   = 2 / abs(  2-phi-sqrt( phi^2-4*phi )  );
    parameters.TOL   = 1;    
    parameters.MITE  = 60;
    parameters.MSAT  = parameters.MITE / 20;
    parameters.SEED  = -1;
    
    parameters.knownMin = -inf;
    parameters.numFeatures = 3;
    parameters.numRules = 2; % 3 features, 2 rules
    parameters.Nrep = 1;
    parameters.printMode = true;
end

printMode = parameters.printMode;

Nrep        = parameters.Nrep;
Nd          = (parameters.numFeatures + 1) * parameters.numRules;
numRules    = parameters.numRules;
fObjName    = func2str(fObj);
%  fObjBnd     = reshape( [zeros(1,Nd) ones(1,Nd)],[],2 );
fObjBnd     = reshape( [zeros(1,Nd-numRules) zeros(1,numRules) ones(1,Nd-numRules) ones(1,numRules)],[],2 );

%% Run with config
t = tic;
% --- Variable parameters
for n = 1 : Nrep    
    [Pg,fPg,details] = UPSO(fObj , fObjBnd, parameters);
%     [Pg,fPg,details] = ParUPSO(fObj , fObjBnd, parameters);
%     [Pg,fPg,details] = ParUPSOHist(fObj , fObjBnd, parameters);
    %    fprintf("Fitness: %.4g\tTime: %.4g", fPg, details.time);
    %    fprintf("\n");    
end
time = toc(t);
if printMode, fprintf('Time = %.4f\n', time); end
%aa = reshape(Pg,[],numFeatures+1); aa(:,[1 3:4]) *= sum(items); aa(:,2) *= length(items); aa(:,end) = round(aa(:,end)*(numActions-1))
%  aa = reshape(Pg,[],numFeatures+1); aa(:,[1:2]) *= sum(items); aa(:,end) = round(aa(:,end)*(numActions-1));
% Rules = Pg;
end