function [indI,X,Y,U,V,LENG_INDI,ind_tri,ind_quad,i_tri,j_tri,i_quad,j_quad] = ...
                                    MOC_2D_steady_irrotational_extend_expansion(geom,params,indIin,LENG_INDIin,...
                                                                                Xin,Yin,Uin,Vin,...
                                                                                ind_triin,ind_quadin,...
                                                                                i_triin,j_triin,i_quadin,j_quadin)
%
% Direct wall point method for the linear part of the nozzle downstream of the circular arc
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
  
  while 1
     indI ++;
     % Wall point
     J=1;
     LENG_INDI(indI)=1;
     % Point 2: left -running C+    ---> ( J+1 , indI-1 )
     [X(J,indI),Y(J,indI),U(J,indI),V(J,indI)] = MOC_2D_steady_irrotational_wall ( X(J+1,indI-1),Y(J+1,indI-1),U(J+1,indI-1),V(J+1,indI-1),...
                                                                                   geom,params) ;
     if (X(J,indI) > geom.xe) 
       break;
     endif
     ind_tri ++;
     i_tri(1:3,ind_tri) = [ indI-1   ; indI-1 ; indI ] ;
     j_tri(1:3,ind_tri) = [ J    ; J+1    ; J    ] ;
     
     J=2;
     for JJ = 2:LENG_INDI(indI-1)-1 % Internal point
       intersection=false;
     % Point 1: right-running C-    ---> ( J-1 , indI   )
     % Point 2: left -running C+    ---> ( J+1 , indI-1 )
       try
         [xtmp,ytmp,utmp,vtmp] = MOC_2D_steady_irrotational_internal_point( X(JJ+1,indI-1),Y(JJ+1,indI-1),U(JJ+1,indI-1),V(JJ+1,indI-1),...
                                                                            X(J-1,indI  ),Y(J-1,indI  ),U(J-1,indI  ),V(J-1,indI  ),...
                                                                            geom,params) ;
       catch
          % Possible intersection of characteristics -> shock -> solution not supported by the method of characteristics
          % Stop the C- characteristic and copy the values from the upstream C- characteristic that intersects
          xtmp = X(JJ+1,indI-1);        ytmp = Y(JJ+1,indI-1);
          utmp = U(JJ+1,indI-1);        vtmp = V(JJ+1,indI-1);
       end
       
       if ( xtmp < X(J-1,indI) )
         disp('Intersection of left-running characteristics inside the diverging region -> weak shock ')
         ind_tri ++;
         i_tri(1:3,ind_tri) = [ indI ; indI-1 ; indI-1 ] ;
         j_tri(1:3,ind_tri) = [ J-1  ; JJ     ; J+1    ] ;
         continue; % Delete the current point and continue with the next one
       endif
       X(J,indI)=xtmp; Y(J,indI)=ytmp; U(J,indI)=utmp; V(J,indI)=vtmp;
       
       test1 = (X(J,indI)-X(J-1,indI))*(Y(J  ,indI-1)-Y(J-1,indI)) - (Y(J,indI)-Y(J-1,indI))*(X(J  ,indI-1)-X(J-1,indI)) ;
       test2 = (X(J,indI)-X(J-1,indI))*(Y(J+1,indI-1)-Y(J-1,indI)) - (Y(J,indI)-Y(J-1,indI))*(X(J+1,indI-1)-X(J-1,indI)) ;
       if (test1*test2<=0) % Intersection of characteristics
          intersection = true;
          break;
       endif
       
       LENG_INDI(indI) ++;
       ind_quad ++;
       i_quad(1:4,ind_quad) = [ indI   ; indI-1 ; indI-1 ; indI ] ;
       j_quad(1:4,ind_quad) = [ J-1    ; JJ     ; JJ+1   ; J    ] ;
       J++;
     endfor
     
     if (intersection)
      % Copy the rest of the data for the C- characteristic from the upstream C- characteristic
        loclen = LENG_INDI(indI-1) - J ;
        X(J:J+loclen-1,indI) = X(J+1:LENG_INDI(indI-1),indI-1);
        Y(J:J+loclen-1,indI) = Y(J+1:LENG_INDI(indI-1),indI-1);
        U(J:J+loclen-1,indI) = U(J+1:LENG_INDI(indI-1),indI-1);
        V(J:J+loclen-1,indI) = V(J+1:LENG_INDI(indI-1),indI-1);
        LENG_INDI(indI) = LENG_INDI(indI) + loclen ;
        
        disp('Intersection of right-running characteristics inside the diverging region -> weak shock')
        
        ind_tri ++;
        i_tri(1:3,ind_tri) = [ indI   ; indI-1 ; indI ] ;
        j_tri(1:3,ind_tri) = [ J-1    ; J    ; J    ] ;
     else
     % Point on the axis of symmetry
        [X(J,indI),Y(J,indI),U(J,indI),V(J,indI)] = MOC_2D_steady_irrotational_internal_point( X(J-1,indI),-Y(J-1,indI),U(J-1,indI),-V(J-1,indI),...
                                                                                               X(J-1,indI), Y(J-1,indI),U(J-1,indI), V(J-1,indI),...
                                                                                               geom,params) ;
        Y(J,indI) += 1.e-6 ; % To avoid singularity on axis
        LENG_INDI(indI)++;
        ind_tri ++;
        i_tri(1:3,ind_tri) = [ indI   ; indI-1 ; indI ] ;
        j_tri(1:3,ind_tri) = [ J-1    ; LENG_INDI(indI-1)    ; J    ] ;
     endif
  end

  % The last point computed in the last loop is located downstream of the nozzle exit lip point
  % Compute the right-running C- from the nozzle exit lip point with the help of the inverse wall method
  [abc,y4,theta4] = MOC_2D_steady_irrotational_get_geometry(geom.xe,2,geom) ;
  J = 1;
  theta4 = atand(theta4);
  xtmp(1:4,1) = [ X(J+1,indI-1) ; 0. ; X(J,indI-1) ; geom.xe ];
  ytmp(1:4,1) = [ Y(J+1,indI-1) ; 0. ; Y(J,indI-1) ; y4 ];
  utmp(1:4,1) = [ U(J+1,indI-1) ; 0. ; U(J,indI-1) ; 0.                ];
  vtmp(1:4,1) = [ V(J+1,indI-1) ; 0. ; V(J,indI-1) ; 0.                ];
  [U(J,indI),V(J,indI),~,~,~,~] = MOC_2D_steady_irrotational_wall_inverse ( xtmp,ytmp,utmp,vtmp,...
                                                                            theta4,geom,params) ;
  X(J,indI) = geom.xe;
  Y(J,indI) = y4;

  for J = 2:LENG_INDI(indI-1)
     % Point 1: right-running C-    ---> ( J-1 , indI   )
     % Point 2: left -running C+    ---> ( J   , indI-1 )
     [X(J,indI),Y(J,indI),U(J,indI),V(J,indI)] = MOC_2D_steady_irrotational_internal_point( X(J  ,indI-1),Y(J  ,indI-1),U(J  ,indI-1),V(J  ,indI-1),...
                                                                                            X(J-1,indI  ),Y(J-1,indI  ),U(J-1,indI  ),V(J-1,indI  ),...
                                                                                            geom,params) ;
     LENG_INDI(indI) ++;
     ind_quad ++;
     i_quad(1:4,ind_quad) = [ indI   ; indI-1 ; indI-1 ; indI ] ;
     j_quad(1:4,ind_quad) = [ J-1    ; J-1    ; J      ; J    ] ;
  endfor
  % Point on axis of symmetry
  J++;
  LENG_INDI(indI)++;
  [X(J,indI),Y(J,indI),U(J,indI),V(J,indI)] = MOC_2D_steady_irrotational_internal_point( X(J-1,indI),-Y(J-1,indI),U(J-1,indI),-V(J-1,indI),...
                                                                                         X(J-1,indI), Y(J-1,indI),U(J-1,indI), V(J-1,indI),...
                                                                                         geom,params) ;
  Y(J,indI) += 1.e-6 ; % To avoid singularity on axis
  ind_tri ++;
  i_tri(1:3,ind_tri) = [ indI   ; indI-1             ; indI ] ;
  j_tri(1:3,ind_tri) = [ J-1    ; LENG_INDI(indI-1)  ; J    ] ;
endfunction