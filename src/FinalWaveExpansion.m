function deltaP = FinalWaveExpansion(deflection,gamma,nu1,Pstag,pambiant)
  
  [M2] = get_Mach_from_nu(nu1+deflection,gamma) ;
  
  p2static = Pstag / ( 1. + 0.5*(gamma-1)*M2.^2 )^(gamma/(gamma-1)) ;
  
  deltaP = pambiant - p2static ;
  
endfunction
