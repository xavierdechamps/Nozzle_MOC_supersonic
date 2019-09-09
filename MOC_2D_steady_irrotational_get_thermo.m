function [Mach,pressure,density,temperature,sound] = MOC_2D_steady_irrotational_get_thermo(u,v,params)
% Get the thermodynamic data from the velocity components
  CP          = params.gamma * params.R / ( params.gamma - 1 );
  Q           = sqrt( u.^2 + v.^2 ) ;
  temperature = params.T - 0.5 * Q.^2 / CP ;
  sound       = sqrt ( params.gamma * params.R * temperature ); 
  Mach        = Q ./ sound ;
  pressure    = params.P .* ( temperature / params.T ).^(params.gamma/(params.gamma-1)) ;
  density     = pressure ./ ( params.R * temperature ) ;
endfunction
