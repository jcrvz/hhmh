function Y = Rastrigin(X)

% X must be a vector
X = X(:);

% Objective function definition
Nd = size(X);
term1 = X.^2;
term2 = 10*cos(2*pi*X);

f_x = 0;
f_x = f_x + sum(term1 - term2);
f_x = f_x + 10 * Nd;

Y = f_x;
end