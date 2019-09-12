function [xo,yo,uo,vo] = MOC_2D_steady_irrotational_internal_point ( xp,yp,up,vp,...
                                                                     xm,ym,um,vm,...
                                                                     params )
% This function computes the intersection of a left-running C+ characteristic
% and a right-running C- characteristic.
% xp,yp,up,vp are conditions from the left -running characteristic C+
% xm,ym,um,vm are conditions from the right-running characteristic C-
   x_orig = [ xp , xm ] ; y_orig = [ yp , ym ] ;
   u_orig = [ up , um ] ; v_orig = [ vp , vm ] ;
   
   step_current = 0;
   step_max = 50;
   eps_pos = 1.e-8;
   eps_vel = 1.e-5;

% Predictor step
% Used to provide a first estimation of the intersection of the 2 characteristics
   [xo,yo,uo,vo] = MOC_2D_steady_irrotational_solve_internal_point ( xp,yp,up,vp,...
                                                                     xm,ym,um,vm,...
                                                                     x_orig,y_orig,u_orig,v_orig,...
                                                                     params ) ;
   while 1
      step_current++;
% Corrector step
      xcp = xp; ycp = 0.5*(yp+yo) ; ucp = 0.5*(up+uo) ; vcp = 0.5*(vp+vo) ;
      xcm = xm; ycm = 0.5*(ym+yo) ; ucm = 0.5*(um+uo) ; vcm = 0.5*(vm+vo) ;
     [xn,yn,un,vn] = MOC_2D_steady_irrotational_solve_internal_point ( xcp,ycp,ucp,vcp,...
                                                                       xcm,ycm,ucm,vcm,...
                                                                       x_orig,y_orig,u_orig,v_orig,...
                                                                       params ) ;
      error_pos = max( [ xo-xn , yo-yn ] );
      error_vel = max( [ uo-un , vo-vn ] );
      xo = xn ;      yo = yn ;
      uo = un ;      vo = vn ;
      % Check if we converged on the position and on the velocity components
      if (abs(error_pos)<eps_pos && abs(error_vel)<eps_vel)
        break;
      endif
      if (step_current>step_max)
        error('The maximum of iterations for the predictor-corrector algorithm has been reached. Stopping the execution...')
      endif
   end
   
endfunction

function [xn,yn,un,vn] = MOC_2D_steady_irrotational_solve_internal_point ( xp,yp,up,vp,...
                                                                           xm,ym,um,vm,...
                                                                           x_orig,y_orig,u_orig,v_orig,...
                                                                           params )
  
  xo = [xp xm] ;
  yo = [yp ym] ;
  uo = [up um] ;
  vo = [vp vm] ;
  
  V     = sqrt  ( uo.^2 + vo.^2 ) ;
  a     = get_speed_sound( params, V ) ;
  theta = atand (vo./uo) ;
  alpha = asind ( a./V ) ;
  lambda(:,1) = tand ( theta(1) + alpha(1) );  % lambda+
  lambda(:,2) = tand ( theta(2) - alpha(2) );  % lambda-
  Q     = uo.^2 - a.^2 ;
  R     = 2*uo.*vo - Q.*lambda ;
  S     = a.^2 .* vo ./ yo ;
  
% Solve the system matrix to get the position of the intersection of the C+ and C- characteristics
  matA    = [ -lambda(1) , 1 ; ...
              -lambda(2) , 1 ] ;
  vecB    = [ y_orig(1)-lambda(1)*x_orig(1) ; ...
              y_orig(2)-lambda(2)*x_orig(2)  ] ;
  new_pos = matA\vecB ;
  xn   = new_pos(1);
  yn   = new_pos(2);
        
  T       = S(:).*(xn-x_orig(:)) + Q(:).*u_orig(:) + R(:).*v_orig(:) ;
  matA    = [ Q(1) , R(1) ; Q(2) , R(2) ] ;
  new_vel = matA\T ;
  un   = new_vel(1) ;
  vn   = new_vel(2) ;
  
endfunction

function a = get_speed_sound(params,V)
% Speed of sound a² = a0²       - 0.5*(gamma-1)*V²
%                   = gamma*R*T - 0.5*(gamma-1)*V²
  a = sqrt ( params.gamma * params.R * params.T - ...
             0.5 * ( params.gamma - 1.) * V.^2 ) ;
endfunction
