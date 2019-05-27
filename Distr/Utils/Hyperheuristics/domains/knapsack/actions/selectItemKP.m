% function [item, itemID] = knapsackStrategy(mode, items)
function [itemID] = selectItemKP(items, mode)
  switch mode
    case 1
      itemID = 1;      % Default
    case 2
      [~, itemID] = max(items(:,2)); % MaP
    case 4 % Original is 3 - Pending update (testing)
      [~, itemID] = max(items(:,2)./items(:,1)); % MPW
    case 3 % Original is 4 - Pending update (testing)
      [~, itemID] = min(items(:,1)); % MiW
  end
%   item = items(itemID,:);
end  