function [nu,mu] = get_prandtl_meyer_function(Mach,gamma)
% Returns the Prandtl-Meyer function nu and the Mach angle mu
  
  assert(Mach>0); 
  mu = asind(1./Mach);
  
  sqroot = sqrt((gamma+1)/(gamma-1));
  nu = sqroot * atand( sqrt(Mach.^2-1)/sqroot ) - atand( sqrt(Mach.^2-1) );
  
endfunction
