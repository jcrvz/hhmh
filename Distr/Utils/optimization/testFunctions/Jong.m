function Y = Jong(X)  
% X must be a vector
X = X(:);

  Y = sum(X.^2); % Sweeps columns and sum
end