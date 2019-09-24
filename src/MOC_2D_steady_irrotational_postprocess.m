function [fig1,fig2,fig3] = MOC_2D_steady_irrotational_postprocess(geom,params,plots,indI,LENG_INDI,...
                                                                   X,Y,U,V,...
                                                                   ind_tri,ind_quad,...
                                                                   i_tri,j_tri,i_quad,j_quad,LOC_PATCHES)
%
% Post-process the results
%
  disp('Postprocessing the results...')
  
  [Mach,pressure,density,temperature,sound] = MOC_2D_steady_irrotational_get_thermo(U,V,params);

  fig1=figure(1);
  plot(X(1,1:indI)/geom.yt,Y(1,1:indI)/geom.yt,'k','linewidth',2); grid on; hold on;
  xlabel('Axial location x / y_t'); ylabel('Radial location y / y_t')

  if (plots.patches>0)
     plot_colours = true;
     switch plots.patches_data
       case 0
          plot_colours = false;
          output_patch = Mach;
       case 1
          output_patch = Mach;
       case 2
          output_patch = pressure;
       case 3
          output_patch = temperature;
       case 4
          output_patch = density;
       case 5
          output_patch = sqrt(U.^2 + V.^2);
       case 6
          output_patch = atand(V./U);
       otherwise
          error('Unknown type of data to display')
     endswitch
     
     indP = 0; % to update patch_x_tri, patch_y_tri, patch_c_tri
     indP2 = 0; % to locate in LOC_PATCHES
     indN1 = 1; % to locate the nodes
     indN2 = 1; % to locate the nodes
     indF1 = 1; % to locate the faces
     indF2 = 1; % to locate the faces
     colours = [0 0 0;1 0 0;0 0 1;1 0 1];
     for I=1:ind_tri
       nodes_tri(indN1:indN1+2,1:2) = [ X(j_tri(1,I),i_tri(1,I))  ,  Y(j_tri(1,I),i_tri(1,I)) ;...
                                        X(j_tri(2,I),i_tri(2,I))  ,  Y(j_tri(2,I),i_tri(2,I)) ;...
                                        X(j_tri(3,I),i_tri(3,I))  ,  Y(j_tri(3,I),i_tri(3,I))   ] ; % [x y]
       if ( I == LOC_PATCHES(indP2+1,1)+1 )
         indP2++; % change of region [ initial-value line ; initial expansion ; nozzle ; plume ]
       endif
       if ( sum( nodes_tri(indN1:indN1+2,1)<plots.patches_xlim ) == 3 )
         % Plot only patches with abscissae below the maximum plot_patches_xlim chosen by user
         faces_tri(indF1,1:3) = [ indN1:indN1+2 ];
         facesC_tri(indF1,1:3) = colours(indP2+1,:);
         nodesC_tri(indN1:indN1+2,1) = [ output_patch(j_tri(1,I),i_tri(1,I)) ; ...
                                         output_patch(j_tri(2,I),i_tri(2,I)) ; ...
                                         output_patch(j_tri(3,I),i_tri(3,I))   ] ;
         indF1 ++;
       endif
       indN1 += 3;
     endfor
     
     indP = 0;
     indP2 = 0;
     for I=1:ind_quad
       nodes_quad(indN2:indN2+3,1:2) = [ X(j_quad(1,I),i_quad(1,I))  ,  Y(j_quad(1,I),i_quad(1,I)) ;...
                                         X(j_quad(2,I),i_quad(2,I))  ,  Y(j_quad(2,I),i_quad(2,I)) ;...
                                         X(j_quad(3,I),i_quad(3,I))  ,  Y(j_quad(3,I),i_quad(3,I)) ;...
                                         X(j_quad(4,I),i_quad(4,I))  ,  Y(j_quad(4,I),i_quad(4,I))   ] ; % [x y]
       if ( I == LOC_PATCHES(indP2+1,2)+1 )
         indP2++; % change of region [ initial-value line ; initial expansion ; nozzle ; plume ]
       endif
       if ( sum( nodes_quad(indN2:indN2+3,1)<plots.patches_xlim ) == 4 )
         % Plot only patches with abscissae below the maximum plot_patches_xlim chosen by user
         faces_quad(indF2,1:4) = [ indN2:indN2+3 ];
         facesC_quad(indF2,1:3) = colours(indP2+1,:);
         nodesC_quad(indN2:indN2+3,1) = [ output_patch(j_quad(1,I),i_quad(1,I)) ; ...
                                          output_patch(j_quad(2,I),i_quad(2,I)) ; ...
                                          output_patch(j_quad(3,I),i_quad(3,I)) ; ...
                                          output_patch(j_quad(4,I),i_quad(4,I))   ] ;
         indF2 ++;
       endif
       indN2 += 4;
     endfor
     
     if (plot_colours)
     % Plot the left- and right-running characteristics with data
##       patch('Faces',faces_tri,'Vertices',nodes_tri,'FaceVertexCData',nodesC_tri,...
##             'FaceColor','interp','LineWidth',1,'LineStyle','-');
##       patch('Faces',faces_quad,'Vertices',nodes_quad,'FaceVertexCData',nodesC_quad,...
##             'FaceColor','interp','LineWidth',1,'LineStyle','-');
     else
     % Plot the left- and right-running characteristics without data
       patch('Faces',faces_tri,'Vertices',nodes_tri/geom.yt,'FaceVertexCData',facesC_tri,...
             'EdgeColor','flat','FaceColor','none','LineWidth',1,'LineStyle','-');
       patch('Faces',faces_quad,'Vertices',nodes_quad/geom.yt,'FaceVertexCData',facesC_quad,...
             'EdgeColor','flat','FaceColor','none','LineWidth',1,'LineStyle','-');
     endif
     colormap(jet(256)) ;
     colorbar; colorlim = caxis;
     
     % Do not plot the left- and right-running characteristics
       patch('Faces',faces_tri,'Vertices',nodes_tri/geom.yt,'FaceVertexCData',nodesC_tri,...
             'FaceColor','interp','LineStyle','none');
       patch('Faces',faces_quad,'Vertices',nodes_quad/geom.yt,'FaceVertexCData',nodesC_quad,...
             'FaceColor','interp','LineStyle','none');

  endif

  axis equal;
  %axis([0 geom.xe 0 max(max(Y(1,:)))])
  %% Plot some arrows representing the flow direction at each point of intersection
  %quiver(X,Y,U,V)
  %saveas(fig1,'Characteristics.pdf')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot the static pressure on the wall + on the axis + comparison with 
  % 1D theoretical curve when the nozzle is chocked
  fig2=figure(2);
  condplot=zeros(size(X,1),size(X,2)); condplot = logical(condplot);
  wallnodesJ = geom.NI+1 : size(X,2) ;
  wallnodesJ = wallnodesJ ( X(1,wallnodesJ) < geom.xe ) ;
  for I=1:size(X,2) % Find the last index for each column -> this is the axis point
    axisnodesJ = find(U(:,I),1,'last');
    condplot( axisnodesJ , I ) = true;
  endfor
  
  [chocked] = get_Laval_theory(geom,params,X(1,wallnodesJ),Y(1,wallnodesJ)) ;
  semilogy(X(1,wallnodesJ)/geom.yt,pressure(1,wallnodesJ)/params.P,'r','linewidth',2); hold on; % Wall nodes
  semilogy(X(condplot)/geom.yt,pressure(condplot)/params.P,'b','linewidth',2); % Axis nodes
  semilogy(X(1,wallnodesJ)/geom.yt,chocked.pressure,'k','linewidth',2); % 1D theory on axis
  grid on; ylabel('Static pressure / Stagnation pressure [-]'); xlabel('Axial location x / y_t');
  legend('Wall','Axis','1D');
  axis([0 geom.xe 0.01 1]);
  set(gca,'xtick',[0:1:10]);
  fig = gcf;
  set(fig,'PaperUnits','normalized');
  set(fig,'PaperPosition',[0 0 1 0.4]);
  %saveas(fig2,'Pressure.pdf')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot the Mach number on the wall + on the axis + comparison with 
  % 1D theoretical curve when the nozzle is chocked
  fig3=figure(3);
  plot(X(1,wallnodesJ)/geom.yt,Mach(1,wallnodesJ),'r','linewidth',2); hold on; % Wall nodes
  plot(X(condplot)/geom.yt,Mach(condplot),'b','linewidth',2); % Axis nodes
  plot(X(1,wallnodesJ)/geom.yt,chocked.mach,'k','linewidth',2); % 1D theory on axis
  grid on; ylabel('Mach number [-]'); xlabel('Axial location x / y_t');
  legend('Wall','Axis','1D');
  axis([0 geom.xe 1 4]);
  set(gca,'xtick',[0:1:10]);
  fig = gcf;
  set(fig,'PaperUnits','normalized');
  set(fig,'PaperPosition',[0 0 1 0.4]);
  lgd=legend; set(legend,'Location','southeast');
  %saveas(fig3,'Mach_number.pdf')

endfunction