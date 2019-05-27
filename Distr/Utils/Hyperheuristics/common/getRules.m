function [Rules, Actions] = getRules(X, nbRules, nbActions)    
% Rules are obtained by reshaping the design vector defining the nbRules as 
% the % number of rows
% Actions are calculated by 

  Rules   = reshape(X, nbRules, []);
%  Actions = round(Rules(:,end) * (nbActions-1));

  % Assigns actions to [0,nbActions]
  Actions_ = floor(Rules(:,end) * (nbActions)); 
  
  % Fix the case when rand==1
  Actions_(Actions_ == nbActions) = nbActions - 1; 
  
  % Fixes range to [1, nbActions]
  Actions = Actions_ + 1; 
  Rules   = Rules(:,1:end-1);
end