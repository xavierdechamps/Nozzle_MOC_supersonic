function [Mach,pressure,density,temperature,sound] = MOC_2D_steady_irrotational_get_thermo(u,v,params)
% Get the thermodynamic data from the velocity components
##  CP          = params.gamma * params.R / ( params.gamma - 1 );
##  Q           = sqrt( u.^2 + v.^2 ) ;
##  temperature = params.T - 0.5 * Q.^2 / CP ;
##  sound       = sqrt ( params.gamma * params.R * temperature ); 
##  Mach        = Q ./ sound ;
##  pressure    = params.P .* ( temperature / params.T ).^(params.gamma/(params.gamma-1)) ;
##  density     = pressure ./ ( params.R * temperature ) ;
  
  % Square of velocity amplitude
  Q           = u.^2 + v.^2 ;
  
  % Stagnation acoustic speed
  sound_stag  = sqrt(params.gamma * params.R * params.T) ;
  
  % Local acoustic speed
  sound       = sqrt(sound_stag^2 - 0.5*(params.gamma - 1)*Q);
  
  % Mach number
  Mach        = sqrt(Q)./sound;
  
  % Static temperature
  temperature = params.T ./ ( 1 + 0.5*(params.gamma-1)*Mach.^2 ) ;
  
  % Static pressure
  pressure    = params.P ./ ( 1 + 0.5*(params.gamma-1)*Mach.^2 ).^(params.gamma/(params.gamma-1)) ;
  
  % Density
  density     = pressure ./ ( params.R * temperature ) ;
endfunction
