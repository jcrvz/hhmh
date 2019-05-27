function [knapsack] = knapsackFitness(knapsack)
  knapsack.profit = sum(knapsack.items(:,2));
  knapsack.weight = sum(knapsack.items(:,1));
  if knapsack.weight <= knapsack.maxWeight
      knapsack.isValid = true;
  else
      knapsack.isValid = false;
  end
end