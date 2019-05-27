function [knapsack, itemsKP, exFlag] = storeKPItem(knapsack, itemsKP, itemID)
    exFlag = 0;
    if knapsack.freeCapacity < itemsKP(itemID,1)
        exFlag = -1;
        return;
    end
    [knapsack, itemsKP] = updateKP(knapsack, itemsKP, itemID);