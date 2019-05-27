function concentration = GetConcentration(agents, mark, ranges)
  % agents is a matrix. One row per agent. One column per dimension
  % mark is a row vector with info of best agent
  % ranges is a row vector with range per dimension
  Na = size(agents,1);
  Nd = size(agents,2);  
  concentration = 0;
  for idx = 1 : Nd
			ratio = 0;
			for idy = 1 : Na
				ratio = ratio + abs(mark(idx) - agents(idy, idx));
			end
			ratio = ratio / (Na*ranges(idx));
			concentration = concentration + ratio;
	end
	concentration = concentration / Nd;  
end
      