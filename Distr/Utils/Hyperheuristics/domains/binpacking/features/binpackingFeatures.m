% ... WIP ... UNUSED ... Not yet implemented
function [feature_] = binpackingFeatures(mode_, items)
    maxValue = 1;
    if mode_ <= 3
      maxValue = max(items(:,1));    
    elseif mode_ <= 6
      maxValue = max(items(:,2));
    end
    
    if maxValue == 0
        maxValue = 1;
    end
    
    switch mode_
      case 1 % Normalized mean weight
        feature_ = mean(items(:,1)) / maxValue;
      case 2 % Normalized median weight
        feature_ = median(items(:,1)) / maxValue;
      case 3 % Normalized standard deviation (weight)
        feature_ = std(items(:,1)) / maxValue;
      case 4 % Normalized mean profit
        feature_= mean(items(:,2)) / maxValue;
      case 5 % Normalized median profit
        feature_ = median(items(:,2)) / maxValue;
      case 6 % Normalized standard deviation (profit)
        feature_ = std(items(:,2)) / maxValue;
      case 7 % Correlation (weight - profit)
        feature_ = corr(items(:,1), items(:,2)) /2 + 0.5;
    end       
end  