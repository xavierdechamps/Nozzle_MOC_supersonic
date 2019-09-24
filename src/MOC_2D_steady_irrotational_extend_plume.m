function [indI,X,Y,U,V,LENG_INDI,ind_tri,ind_quad,i_tri,j_tri,i_quad,j_quad] = ...
                                    MOC_2D_steady_irrotational_extend_plume(geom,params,indIin,LENG_INDIin,...
                                                                            Xin,Yin,Uin,Vin,...
                                                                            ind_triin,ind_quadin,...
                                                                            i_triin,j_triin,i_quadin,j_quadin)
%
% Check the static pressure at the nozzle lip point. If static pressure > ambiant static pressure: Prandtl-Meyer expansion
%                                                    If static pressure < ambiant static pressure: weak shock
% The Prandtl-Meyer expansion is discretized by geom.NIexpansion lines.
% The free-pressure point method is used to obtain the shape of the plume.
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
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [Mach,pressure,~,~,~] = MOC_2D_steady_irrotational_get_thermo(U,V,params);
  plume_case = 1; % [1/2] 1: weak shock at exit lip point
                  %       2: Prandtl-Meyer expansion at exit lip point
  if ( params.PRatio < 1 )
  %if ( pressure(1,end) < params.Pstatic )
    disp('Oblique shock wave emanating from the exit lip point')
    if ( params.PRatio > 0.5 )
      disp('Weak shock wave -> MOC can continue in the plume without special treatment')
    else
      disp('Strong shock wave at the lip point. The MOC cannot handle such problems.')
      disp('Stopping the computation now...')
      return;
    endif
  else
    disp('Prandtl-Meyer expansion wave at the exit lip point')
    disp('Need to iterate to find the final angle that matches the static pressure outside the nozzle.')
    plume_case = 2;
  endif

  params.Pstatic = pressure(1,indI) / params.PRatio ;

  % At this position of the code, the plume can be calculated with the method of characteristics
  
  if ( plume_case == 2 ) % Prandtl-Meyer expansion
   
   % We need to iterate on the final angle theta so that the static pressure matches the ambiant static pressure
     [nu1,mu1] = get_prandtl_meyer_function(Mach(1,indI),params.gamma); % Mach angles before the expansion
     % Initial values, to launch the iterative process
     deflection = 30; % deflection angle in degrees nu2 = nu1 + deflection
     disp(['    Mach number before the expansion: ',num2str(Mach(1,indI))])
     % Find the deflection that matches the ambiant static pressure and the static pressure after the expansion
     deflection = fzero(@(delta12) FinalWaveExpansion(delta12,params.gamma,nu1,params.P,params.Pstatic),[0 deflection],optimset('TolX',1e-12));
     Mach2      = get_Mach_from_nu(nu1+deflection,params.gamma);
     disp(['    Mach number after the expansion: ',num2str( Mach2 )      ])
     disp(['    Deflection angle (degrees)       ',num2str( deflection ) ])
     
     incrdefl = deflection / geom.NIexpansion ;
     for I=1:geom.NIexpansion-1 % from 1st expansion wave to penultimate expansion wave
      % Get data for exit lip point for current expansion wave
        Mach2 = get_Mach_from_nu(nu1+incrdefl*I,params.gamma);
        [~,~,~,a2] = MOC_2D_steady_irrotational_get_thermo2(Mach2,params) ;
        J = 1;
        indI++;
        LENG_INDI(indI)=1;
        X(J,indI) = X(J,indI-1);  Y(J,indI) = Y(J,indI-1); 
        U(J,indI) = Mach2 * a2 * cosd(geom.te+incrdefl*I) ; 
        V(J,indI) = Mach2 * a2 * sind(geom.te+incrdefl*I) ;
        
       % Get data for the internal points
        for J = 2:LENG_INDI(indI-1)
        % Point 1: right-running C-    ---> ( J-1 , I   )
        % Point 2: left -running C+    ---> ( J   , I-1 )
           [X(J,indI),Y(J,indI),U(J,indI),V(J,indI)] = MOC_2D_steady_irrotational_internal_point( X(J  ,indI-1),Y(J  ,indI-1),U(J  ,indI-1),V(J  ,indI-1),...
                                                                                                  X(J-1,indI  ),Y(J-1,indI  ),U(J-1,indI  ),V(J-1,indI  ),...
                                                                                                  geom,params) ;
           LENG_INDI(indI)++;
           ind_quad ++;
           i_quad(1:4,ind_quad) = [ indI  ; indI-1 ; indI-1 ; indI ] ;
           j_quad(1:4,ind_quad) = [ J-1   ; J-1    ; J      ; J ] ;
        endfor
       % Point on the axis of symmetry
        J++;
        [X(J,indI),Y(J,indI),U(J,indI),V(J,indI)] = MOC_2D_steady_irrotational_internal_point( X(J-1,indI),-Y(J-1,indI),U(J-1,indI),-V(J-1,indI),...
                                                                                               X(J-1,indI), Y(J-1,indI),U(J-1,indI), V(J-1,indI),...
                                                                                               geom,params) ;
        Y(J,indI) += 1.e-6 ; % To avoid singularity on axis
        LENG_INDI(indI)++;
        ind_tri ++;
        i_tri(1:3,ind_tri) = [ indI ; indI-1            ; indI ] ;
        j_tri(1:3,ind_tri) = [ J-1  ; LENG_INDI(indI-1) ; J ] ;
        
     endfor
  endif
  
 % At this point, we are either at the final stage of the Prandtl-Meyer expansion or just after the weak shock
  while 1
      indI ++;
      % Jet point
      J=1;
      LENG_INDI(indI)=1;
      % Point 2: left -running C+    ---> ( J+1 , indI-1 )indI ++;
      [X(J,indI),Y(J,indI),U(J,indI),V(J,indI)] = MOC_2D_steady_irrotational_free_pressure (X(J  ,indI-1),Y(J  ,indI-1),U(J  ,indI-1),V(J  ,indI-1),...
                                                                                            X(J+1,indI-1),Y(J+1,indI-1),U(J+1,indI-1),V(J+1,indI-1),...
                                                                                            geom,params ) ;
      if ( X(J,indI) > (geom.xe+geom.xplume) )
         indI --;
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
          disp('Intersection of left-running characteristics inside the plume -> weak shock ')
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
         
         disp('Intersection of right-running characteristics inside the plume -> weak shock')
         
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
     

endfunction