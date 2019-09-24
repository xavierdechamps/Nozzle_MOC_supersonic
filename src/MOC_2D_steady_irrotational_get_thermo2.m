function [pressure,density,temperature,sound] = MOC_2D_steady_irrotational_get_thermo2(Mach,params)
% Get the thermodynamic data from the Mach number
##  CP          = params.gamma * params.R / ( params.gamma - 1 );
##  Q           = sqrt( u.^2 + v.^2 ) ;
##  temperature = params.T - 0.5 * Q.^2 / CP ;
##  sound       = sqrt ( params.gamma * params.R * temperature );
##  pressure    = params.P .* ( temperature / params.T ).^(params.gamma/(params.gamma-1)) ;
##  density     = pressure ./ ( params.R * temperature ) ;
  
  % Static temperature
  temperature = params.T ./ ( 1 + 0.5*(params.gamma-1)*Mach.^2 ) ;
  
  % Local acoustic speed
  sound       = sqrt( params.gamma * params.R * temperature );  
  
  % Static pressure
  pressure    = params.P ./ ( 1 + 0.5*(params.gamma-1)*Mach.^2 ).^(params.gamma/(params.gamma-1)) ;
  
  % Density
  density     = pressure ./ ( params.R * temperature ) ;
endfunction
