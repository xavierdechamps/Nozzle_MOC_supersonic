function [xsonic,xvnull,u] = MOC_2D_steady_irrotational_IVLINE ( geom , params , y )
  
  % Origin of the coordinate system relative to the nozzle throat [m]
  eps = -geom.yt * sqrt((params.gamma+1)*(geom.delta+1)*geom.yt/geom.rhou) * 0.5 / (3+geom.delta) ;
  
  % Sonic speed a*
  astar = sqrt(2*params.gamma*params.R*params.T/(params.gamma+1)) ;
  
  % Coefficient of the linear axial perturbation velocity [1/m]
  alpha = sqrt((1+geom.delta)/(geom.yt*geom.rhou*(params.gamma+1)));
  
  % Equation for the sonic line x = coeff_sonic * y^2, where coeff_sonic is [1/m]
  coeff = -(params.gamma+1)*alpha*0.5/(1+geom.delta);
  xsonic=coeff*y.^2 - eps;
  
  % Equation for the v=0 line: x=coeff_vnull * y^2, where coeff_vnull is [1/m]
  coeff_vnull = -(params.gamma+1)*alpha*0.5/(3+geom.delta);
  xvnull=coeff_vnull*y.^2 - eps;
  
  % Perturbation velocity field on the v=0 line
  u = astar * (1 + alpha * (xvnull+eps) + (params.gamma+1)*(alpha^2)*(y.^2)*0.5/(1+geom.delta) );
  
endfunction