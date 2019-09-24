function [indI,X,Y,U,V,LENG_INDI,ind_tri,ind_quad,i_tri,j_tri,i_quad,j_quad] = ...
                                    MOC_2D_steady_irrotational_extend_initial_line(geom,params,...
                                                                                   xvnull,ythroat,uvnull,vvnull,...
                                                                                   Xin,Yin,Uin,Vin)
%
% Extend the solution from the initial-value line, fill-in the matrices X, Y, U and V
% which contains the (x,y) coordinates of the intersections of the characteristics
% and the (u,v) velocity components at these points
%
  LENG_INDI(1)   = 1;
  ind_tri  = 0;
  ind_quad = 0;
  X              = Xin;
  Y              = Yin;
  U              = Uin;
  V              = Vin;
  for I = 2:geom.NI
     X(1,I) = xvnull(I); Y(1,I) = ythroat(I); U(1,I) = uvnull(I); V(1,I) = vvnull;
     LENG_INDI(I) = 1;
     for J = 2:2*(I-1) % Internal points
     % Point 1: right-running C-    ---> ( J-1 , I   )
     % Point 2: left -running C+    ---> ( J-1 , I-1 )
        [X(J,I),Y(J,I),U(J,I),V(J,I)] = MOC_2D_steady_irrotational_internal_point( X(J-1,I-1),Y(J-1,I-1),U(J-1,I-1),V(J-1,I-1),...
                                                                                   X(J-1,I  ),Y(J-1,I  ),U(J-1,I  ),V(J-1,I  ),...
                                                                                   geom,params) ;
        LENG_INDI(I)++;
        if (J==2)
          ind_tri ++;
          i_tri(1:3,ind_tri) = [ I   ; I-1 ; I ] ;
          j_tri(1:3,ind_tri) = [ J-1 ; J-1 ; J ] ;
        else
          ind_quad ++;
          i_quad(1:4,ind_quad) = [ I   ; I-1 ; I-1 ; I ] ;
          j_quad(1:4,ind_quad) = [ J-1 ; J-2 ; J-1 ; J ] ;
        endif
     endfor
     % Point on the axis of symmetry
     J++;
     [X(J,I),Y(J,I),U(J,I),V(J,I)] = MOC_2D_steady_irrotational_internal_point( X(J-1,I),-Y(J-1,I),U(J-1,I),-V(J-1,I),...
                                                                                X(J-1,I), Y(J-1,I),U(J-1,I), V(J-1,I),...
                                                                                geom,params) ;
     Y(J,I) += 1.e-6 ; % To avoid singularity on axis
     LENG_INDI(I)++;
     ind_tri ++;
     i_tri(1:3,ind_tri) = [ I   ; I-1 ; I ] ;
     j_tri(1:3,ind_tri) = [ J-1 ; J-2 ; J ] ;
  endfor
  
  indI = I;
  
endfunction
