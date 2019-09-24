function [xo,yo,uo,vo] = MOC_2D_steady_irrotational_free_pressure ( xw,yw,uw,vw,xp,yp,up,vp,geom,params )
%
% This function computes the intersection of a left-running C+ characteristic
% with the plume external jet. Output of the function is the position + data at the intersection on the jet.
% xw,yw,uw,vw = point 3 on jet, previous upwards line
% xp,yp,up,vp = point 2 internal point, previous upwards line
%
   x_orig = [xp  xw  xw ] ; y_orig = [yp  yw  yw ] ;
   u_orig = [up  uw  uw ] ; v_orig = [vp  vw  vw ] ;
   
   step_current = 0;
   step_max = 50;
   eps_pos = 1.e-8;
   eps_vel = 1.e-5;
   
% Predictor step
% Used to provide a first estimation of the intersection of the characteristic with the wall
   [xo,yo,uo,vo] = MOC_2D_steady_irrotational_solve_free_pressure ( xp,yp,up,vp,...
                                                                    x_orig,y_orig,u_orig,v_orig,...
                                                                    geom,params ) ;
   % The first estimation can lead to complex numbers, keep only the real part for the corrector step
   xo=real(xo); yo=real(yo); uo=real(uo); vo=real(vo);
   while 1
      step_current++;
% Corrector step
      xcp = xp;  ycp = 0.5*(yp+yo);  ucp = 0.5*(up+uo);  vcp = 0.5*(vp+vo);
     [xn,yn,un,vn] = MOC_2D_steady_irrotational_solve_free_pressure ( xcp,ycp,ucp,vcp,...
                                                                      x_orig,y_orig,u_orig,v_orig,...
                                                                      geom,params ) ;
      error_pos = max( [ xo-xn , yo-yn ] );
      error_vel = max( [ uo-un , vo-vn ] );
      xo = xn ;      yo = yn ;
      uo = un ;      vo = vn ;
      x_orig(3) = xo; % data for point 4, reinjected in MOC_2D_steady_irrotational_solve_free_pressure
      y_orig(3) = yo; 
      u_orig(3) = uo; 
      v_orig(3) = vo; 
      
      % Check if we converged on the position and on the velocity components
      if (abs(error_pos)<eps_pos && abs(error_vel)<eps_vel)
        break;
      endif
      if (step_current>step_max)
        error('The maximum of iterations for the predictor-corrector algorithm has been reached. Stopping the execution...')
      endif
   end
   
endfunction


function [xn,yn,un,vn] = MOC_2D_steady_irrotational_solve_free_pressure ( xp,yp,up,vp,...
                                                                          x_orig,y_orig,u_orig,v_orig,...
                                                                          geom,params )
  
  % Velocity at solution point, flow is isentropic
  V4    = sqrt(2*params.gamma*params.R*params.T*(1-(params.Pstatic/params.P)^((params.gamma-1)/params.gamma))/(params.gamma-1)); % free pressure
  V     = sqrt  ( up.^2 + vp.^2 ) ;
  a     = get_speed_sound( params, V ) ;
  theta = atand (vp./up) ;
  alpha = asind ( a./V ) ;
  lambda = tand ( theta + alpha );  % lambda+
  Q     = up.^2 - a.^2 ;
  R     = 2*up.*vp - Q.*lambda ;
  S     = geom.delta * a.^2 .* vp ./ yp ;
  lambda0 = (v_orig(2)+v_orig(3))/(u_orig(2)+u_orig(3));  % slope of the jet boundary
%  if (step>0)
    % Better estimation because we already have some data for the point to be calculated
%    lambda0 = (v_orig(1,2)+v_orig(2,1))/(u_orig(1,2)+u_orig(2,1)) ; % slope of the jet boundary
%  else
    % First estimation with data at point 3
%    lambda0 = v_orig(2,1)/u_orig(2,1) ; % slope of the jet boundary
%  endif
  
% Solve the system matrix to get the position of the intersection of the C+ characteristic with the jet
%   Left-running C+ characteristic
%      y4 - lambda  * x4 = y_orig(1) - lambda  * x_orig(1)
%   Jet boundary
%      y4 - lambda0 * x4 = y_orig(2) - lambda0 * x_orig(2)
  matA=[ - lambda , 1 ;
         -lambda0 , 1 ];
  RHS = [y_orig(1) - lambda  * x_orig(1) ; y_orig(2) - lambda0 * x_orig(2)] ;
  sol = matA\RHS;
  xn = sol(1);
  yn = sol(2);
    
  T  = S*(xn-x_orig(1)) + Q*u_orig(1) + R*v_orig(1) ;
  un = (Q*T-R*sqrt(V4^2*(Q^2+R^2)-T^2))/(Q^2+R^2);
  vn = sqrt(V4^2-un^2);
  
endfunction