function printRules(Rules, Actions)
  fprintf("--- Rules:\n");
  for idx = 1 : size(Rules,1)
    fprintf("---\t");
    fprintf("\t%6f", Rules(idx,:));
    fprintf("\t-->\t%d\n",Actions(idx));    
  end  
end