function Y = TransformConcentration(X, mode)
% X = original value 
% mode =
%          0  :  Concentration
%          1  :  Spread
%          2  :  Dispersion
  switch (mode)
    case 0
      a = 20e3;
      Y = 2*( (exp(-a*X)-exp(-a)) ./ (1+exp(-a*X)) );
    case 1
%      mu = 0.2;
%      sigma = 0.05;
%      Y = 1/(sqrt(2*pi)*sigma*8) * exp(-(X-mu).^2./(2*sigma^2));
      mu = 0.01;
      sigma = mu / 2;
      Y = 1/(sqrt(2*pi)*sigma*80) * exp(-(X-mu).^2./(2*sigma^2));
    case 2
%      a = 8;
      a = 4;
      Y = 2./(1+exp(-a*(X-1)));
  end
end