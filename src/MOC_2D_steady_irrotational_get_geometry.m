function [abc,y,tangent] = MOC_2D_steady_irrotational_get_geometry(x,step,geom)
% Get the geometry of a nozzle, see page 128, Gas Dynamics V2
% The divergent region contains a circular arc of radius radtd, which extends
% up to thetaa degrees. At the attachment point A located at an angle thetaa,
% a 2nd order poly is tangent to the circular arc and extends to an axial position
% xe. The angle of the 2nd order poly at xe is given by thetae. 
  
  radtu = geom.rhou ; % Radius of curvature of upstream   circular arc
  radtd = geom.rhod ; % Radius of curvature of downstream circular arc
  yt    = geom.yt   ; % Radius of nozzle throat
  thetaa= geom.ta   ; % Attachment angle between downstream arc and 2nd order poly
  thetae= geom.te   ; % Nozzle exit lip angle
  xe    = geom.xe   ; % Nozzle length
  
  % Coordinates of attachment point between downstream arc and 2nd order poly
  xa    = radtd * sind ( thetaa ) ;
  ya    = yt + ( 1 - cosd ( thetaa ) ) * radtd ;
  
  % 2nd order poly       y = a + bx + cx^2
  %                tangent =     b  + 2cx
% Get the coefficients a, b and c for the 2nd order poly
% 1st line of matA: (x,y) at attachment point A are known
% 2nd line of matA: slope at attachment point A is known = tan ( thetaa )
% 3rd line of matA: slope  at exit lip point is known = tan ( thetae )
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