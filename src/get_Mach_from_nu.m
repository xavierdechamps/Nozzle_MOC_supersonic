function [Mach] = get_Mach_from_nu(nu,gamma)
% Returns the Mach number from the Prandtl-Meyer function nu.
% nu is supposed to be given in degrees.

  sqroot = sqrt((gamma+1)/(gamma-1));
  equation = @(Mach,nu) sqroot * atand( sqrt(Mach.^2-1)/sqroot ) - atand( sqrt(Mach.^2-1) ) - nu ;
  
  guessMach = [1.00001 50 ];
  for k=1:length(nu)
    assert (nu(k)>0);
    assert (nu(k)<=125.);
    ma(k) = fzero(@(Mach) equation(Mach,nu(k)), guessMach);
  endfor
  
  Mach = ma;
endfunction
