function [knapsack, itemsKP] = updateKP(knapsack, itemsKP, itemID)        
    knapsack.items = [knapsack.items; itemsKP(itemID,:)];
    knapsack.freeCapacity = knapsack.freeCapacity - itemsKP(itemID,1);
    itemsKP = [itemsKP(1:itemID-1,:); itemsKP(itemID+1:end,:)];