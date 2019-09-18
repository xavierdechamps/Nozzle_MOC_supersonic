function [xn,yn,un,vn] = MOC_2D_steady_irrotational_solve_free_pressure ( xp,yp,up,vp,...
                                                                          x_orig,y_orig,u_orig,v_orig,...
                                                                          pressure4,step,params )
  
  % Velocity at solution point, flow is isentropic
  V4    = sqrt(2*params.gamma*params.R*params.T*(1-(pressure4/params.P)^((params.gamma-1)/params.gamma))/(params.gamma-1)) ;% free pressure
  V     = sqrt  ( up.^2 + vp.^2 ) ;
  a     = get_speed_sound( params, V ) ;
  theta = atand (vp./up) ;
  alpha = asind ( a./V ) ;
  lambda = tand ( theta + alpha ) ; % lambda+
  Q     = up.^2 - a.^2 ;
  R     = 2*up.*vp - Q.*lambda ;
  S     = a.^2 .* vp ./ yp ;
  if (step>0)
    % Better estimation because we already have some data for the point to be calculated
    lambda0 = (v_orig(1,2)+v_orig(2,1))/(u_orig(1,2)+u_orig(2,1)) ; % slope of the jet boundary
  else
    % First estimation with data at point 3
    lambda0 = v_orig(2,1)/u_orig(2,1) ; % slope of the jet boundary
  endif
  
% Solve the system matrix to get the position of the intersection of the C+ characteristic with the wall
%   Left-running C+ characteristic
%      y4 - lambda  * x4 = y_orig(1) - lambda  * x_orig(1)
%   Jet boundary
%      y4 - lambda0 * x4 = y_orig(2) - lambda0 * x_orig(2)
  matA=[ - lambda , 1 ;
         -lambda0 , 1 ];
  RHS = [y_orig(1,1) - lambda  * x_orig(1,1) ; y_orig(2,1) - lambda0 * x_orig(2,1)] ;
  sol = matA\RHS;
  xn = sol(1);
  yn = sol(2);
    
  T  = S*(xn-x_orig(1,1)) + Q*u_orig(1,1) + R*v_orig(1,1) ;
  un = (Q*T-R*sqrt(V4^2*(Q^2+R^2)-T^2))/(Q^2+R^2);
  vn = sqrt(V4^2-un^2);
  
endfunction