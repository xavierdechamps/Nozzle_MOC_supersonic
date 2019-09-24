function [indI,X,Y,U,V,LENG_INDI,ind_tri,ind_quad,i_tri,j_tri,i_quad,j_quad] = ...
                                    MOC_2D_steady_irrotational_extend_initial_expansion(geom,params,indIin,LENG_INDIin,...
                                                                                        Xin,Yin,Uin,Vin,...
                                                                                        ind_triin,ind_quadin,...
                                                                                        i_triin,j_triin,i_quadin,j_quadin)
%
% Inverse wall point method for the downwards circular arc with prespecified points
% which contains the (x,y) coordinates of the intersections of the characteristics
% and the (u,v) velocity components at these points
%
  indI = indIin;
  X    = Xin ; 
  Y    = Yin ;
  U    = Uin ;
  V    = Vin ;
  LENG_INDI = LENG_INDIin;
  ind_tri   = ind_triin;
  ind_quad  = ind_quadin;
  i_tri     = i_triin;
  j_tri     = j_triin;
  i_quad    = i_quadin;
  j_quad    = j_quadin;
  for I = 1:length(geom.circdownX)
     indI ++;
     theta4 = geom.circdownTheta(I);
     LENG_INDI(indI)=1;
     J=1;
     % Prescribed point on the circular arc -> get the values of U and V
     xtmp(1:4,1) = [ X(J+1,indI-1) ; 0. ; X(J,indI-1) ; geom.circdownX(I) ];
     ytmp(1:4,1) = [ Y(J+1,indI-1) ; 0. ; Y(J,indI-1) ; geom.circdownY(I) ];
     utmp(1:4,1) = [ U(J+1,indI-1) ; 0. ; U(J,indI-1) ; 0.                ];
     vtmp(1:4,1) = [ V(J+1,indI-1) ; 0. ; V(J,indI-1) ; 0.                ];
    [U(J,indI),V(J,indI),newx2(I),newy2(I),newu2(I),newv2(I)] = MOC_2D_steady_irrotational_wall_inverse ( xtmp,ytmp,utmp,vtmp,...
                                                                                                          theta4,geom,params) ;
     X(J,indI) = geom.circdownX(I);
     Y(J,indI) = geom.circdownY(I);
     
     J=2;
     for JJ = 2:LENG_INDI(indI-1) % Internal points
     % Point 1: right-running C-    ---> ( J-1 , indI   )
     % Point 2: left -running C+    ---> ( J   , indI-1 )
        [xtmp,ytmp,utmp,vtmp] = MOC_2D_steady_irrotational_internal_point( X(JJ,indI-1),Y(JJ,indI-1),U(JJ,indI-1),V(JJ,indI-1),...
                                                                           X(J-1,indI  ),Y(J-1,indI  ),U(J-1,indI  ),V(J-1,indI  ),...
                                                                           geom,params) ;
        if (ytmp>=Y(1,indI))
         % If the C+ characteristic exits the geometry of the nozzle, delete the last computed point
           disp('... Initial expansion: deleting a point outside of the nozzle')
           continue;
        endif
        X(J,indI)=xtmp; Y(J,indI)=ytmp; U(J,indI)=utmp; V(J,indI)=vtmp;
        LENG_INDI(indI)++;
        ind_quad ++;
        i_quad(1:4,ind_quad) = [ indI   ; indI-1 ; indI-1 ; indI ] ;
        j_quad(1:4,ind_quad) = [ J-1    ; JJ-1    ; JJ      ; J    ] ;
        J++;
     endfor
     % Point on the axis of symmetry
     [X(J,indI),Y(J,indI),U(J,indI),V(J,indI)] = MOC_2D_steady_irrotational_internal_point( X(J-1,indI),-Y(J-1,indI),U(J-1,indI),-V(J-1,indI),...
                                                                                            X(J-1,indI), Y(J-1,indI),U(J-1,indI), V(J-1,indI),...
                                                                                            geom,params) ;
     Y(J,indI) += 1.e-6 ; % To avoid singularity on axis
     LENG_INDI(indI)++;
     ind_tri ++;
     i_tri(1:3,ind_tri) = [ indI   ; indI-1 ; indI ] ;
     j_tri(1:3,ind_tri) = [ J-1    ; JJ    ; J    ] ;
  endfor
endfunction