function Y = ShiftedJong(X)  
% X must be a vector
X = X(:);

A = 2;
Nd = numel(X);
idx = 1:Nd;
  Y = sum((X-idx*A).^2); % Sweeps columns and sum
end