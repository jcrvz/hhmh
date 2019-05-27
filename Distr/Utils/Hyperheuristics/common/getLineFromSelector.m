function [ml, b] = getLineFromSelector(selector)
m = ( selector(2,2) - selector(1,2) ) / (selector(2,1) - selector(1,1));
ml = -1/m;
x = mean(selector,1);
b = x(2) - ml*x(1);
end