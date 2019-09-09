function [abc,y,tangent] = MOC_2D_steady_irrotational_get_geometry(x,step,geom)
% Get the geometry of a nozzle, see page 128, Gas Dynamics V2
  
  radtu = geom.rhou ; %0.05   ; % Radius of curvature of upstream   circular arc
  radtd = geom.rhod ; %0.0125 ; % Radius of curvature of downstream circular arc
  yt    = geom.yt   ; %0.025  ; % Radius of nozzle throat
  thetaa= geom.ta   ; %35     ; % Attachment angle between downstream arc and 2nd order poly
  thetae= geom.te   ; %10     ; % Nozzle exit lip angle
  xe    = geom.xe   ; %0.25   ; % Nozzle length
  
  % Coordinates of attachment point between downstream arc and 2nd order poly
  xa    = radtd * sind ( thetaa ) ;
  ya    = yt + ( 1 - cosd ( thetaa ) ) * radtd ;
  
  % 2nd order poly       y = a + bx + cx^2
  %                tangent = b + 2cx
  matA = [ 1 , xa , xa^2 ;
           0 , 1  , 2*xa ;
           0 , 1  , 2*xe  ] ;
  RHS = [ ya ;
          tand(thetaa) ;
          tand(thetae) ] ;
  abc = matA\RHS ;
  
  y       = 0.;
  tangent = 0.;
  if (step==2) % get the geometry after the downstream circular arc
    assert ( x >= xa );
%    assert ( x <= xe );
    y       = abc(1) + abc(2) * x + abc(3) * x.^2 ;
    tangent = abc(2) + 2*abc(3) * x ;
  endif
endfunction