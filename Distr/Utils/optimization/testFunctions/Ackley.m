function Y = Ackley(X)
% X must be a vector
X = X(:);

% Objective function definition
Nd = numel(X);
sum1 = sum(X.^2);
sum2 = sum(cos(2*pi*X),2);

Y = 20 + exp(1) - 20 * exp(-0.2*sqrt(sum1/Nd))  -  exp( sum2/Nd );
end