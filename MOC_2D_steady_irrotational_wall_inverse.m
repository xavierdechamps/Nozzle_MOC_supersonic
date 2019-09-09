function [u4,v4,x2,y2,u2,v2] = MOC_2D_steady_irrotational_wall_inverse ( xin,yin,uin,vin,...
                                                                         theta4,params )
   
   step_current = 0;
   step_max = 10;
   eps_pos = 1.e-8;
   eps_vel = 1.e-5;
   
   x=xin;y=yin;u=uin;v=vin;
% Predictor step
   x(1,2) = x(2,1);  y(1,2) = y(2,1); u(1,2) = u(2,1) ; v(1,2) = v(2,1) ;
   [x,y,u,v] = MOC_2D_steady_irrotational_solve_wall_inverse ( x,y,u,v,theta4,params) ;
   while 1
      step_current++;
% Corrector step
      x(1,2) = x(2,1);  y(1,2) = y(2,1); u(1,2) = 0.5*(u(1,2)+u(2,1)) ; v(1,2) = 0.5*(v(1,2)+v(2,1)) ;
     [xn,yn,un,vn] = MOC_2D_steady_irrotational_solve_wall_inverse ( x,y,u,v,theta4,params) ;
      
      error_vel = max( [ u(3,1)-un(3,1) , v(3,1)-vn(3,1) ] );
      x = xn ;
      y = yn ;
      u = un ;
      v = vn ;
      if (abs(error_vel)<eps_vel)
        break;
      endif
      if (step_current>step_max)
        error('The maximum of iterations for the predictor-corrector algorithm has been reached. Stopping the execution...')
      endif
    end
    
    u4 = u(3,1);
    v4 = v(3,1);
    x2 = x(1,2);
    y2 = y(1,2);
    u2 = u(1,2);
    v2 = v(1,2);
endfunction

function [xn,yn,un,vn] = MOC_2D_steady_irrotational_solve_wall_inverse ( x,y,u,v,theta4,params )
% x contains the absissae of the 3 nodes required to compute the internal point 2
%    1st node internal (known position and data)
%    2nd node internal (unknown position, unknown data)
%    3rd node at wall (known position and data)
%    4th node at wall (known position, unknown data)
%
  xn=x;yn=y;un=u;vn=v;
%%**** First step, get the position [xp,yp] of the unknown internal point
%   Predictor step, assume that u2,v2 = u3,v3
  eps_vel = 1.e-5;
  steps_max = 10;
  step_current = 1;
  while 1
    V     = sqrt  ( un(1,2)^2 + vn(1,2)^2 ) ;
    a     = get_speed_sound( params, V ) ;
    theta = atand (vn(1,2)/un(1,2)) ;
    alpha = asind ( a/V ) ;
    lambda = tand ( theta + alpha ) ; % lambda+
    % Slope of line 1-3
    slope13 = (y(2,1)-y(1,1))/(x(2,1)-x(1,1));
    matA = [ -slope13 1 ;
             -lambda  1 ];
    RHS  = [ y(1,1)-slope13*x(1,1) ;
             y(3,1)-lambda*x(3,1)  ];
    sol  = matA \ RHS;
    xn(1,2)= sol(1);
    yn(1,2)= sol(2);
    % Interpolation of data at points 1 and 3 to get data at point 2
    unn  = u(1,1)+(xn(1,2)-x(1,1))/(x(2,1)-x(1,1))*(u(2,1)-u(1,1));
    vnn  = v(1,1)+(xn(1,2)-x(1,1))/(x(2,1)-x(1,1))*(v(2,1)-v(1,1));
    if (abs(unn-un(1,2))<eps_vel && abs(vnn-vn(1,2))<eps_vel)
      break;
    endif
    if (step_current>steps_max)
      error('Could not converge to a position for the new internal point')
    endif
    un(1,2) = unn;
    vn(1,2) = vnn;
    step_current++;
  end
  
  % Get the data at point 4 on the wall
  Q     = un(1,2)^2 - a^2 ;
  R     = 2*un(1,2)*vn(1,2) - Q*lambda ;
  S     = a^2 * vn(1,2) / yn(1,2) ;
  T     = S*(x(3,1)-xn(1,2)) + Q*un(1,2) + R*vn(1,2) ;
  matA  = [ tand(theta4) , -1 ;
            Q            ,  R  ];
  RHS   = [ 0 ; T ];
  sol   = matA \ RHS ;
  un(3,1) = sol(1) ;
  vn(3,1) = sol(2) ;
  
endfunction

function a = get_speed_sound(params,V)
% Speed of sound a² = a0²       - 0.5*(gamma-1)*V²
%                   = gamma*R*T - 0.5*(gamma-1)*V²
  a = sqrt ( params.gamma * params.R * params.T - ...
             0.5 * ( params.gamma - 1.) * V^2 ) ;
endfunction
