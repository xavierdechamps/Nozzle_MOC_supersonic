function [chocked] = get_Laval_theory(geom,params,X,Y)
% Get the isentropic evolution of the flow through a Laval nozzle, after the throat section

   gam = params.gamma;
   gamm1=gam-1;
   gamp1=gam+1;

   area = Y.^2 / (Y(1)^2);

   ##xline = geom.rhod*sind( geom.ta )+0.01 : 0.01 : geom.xe  ;
   ##xnozzle = [
   ##      geom.rhod*sind(geom.circdownTheta)  , ...
   ##      xline
   ##      ]   ;
   ##ynozzle = geom.yt + [ 
   ##                geom.rhod*(1-cosd(geom.circdownTheta))  , ...
   ##                geom.rhod*(1-cosd(geom.ta)) + (xline-geom.rhod*sind(geom.ta))*tand(geom.te)  
   ##            ]   ;
   ##area=ynozzle.^2 / (ynozzle(1)^2); % Normalized by area at throat

% Equation linking the area ratio and the Mach number ratio
% m = M/M* and a = A/A* where the star denotes conditions at the throat
   the_eq = @(m,a) ( 2*(1+0.5*gamm1*m^2)/gamp1 )^(0.5*gamp1/gamm1) /m - a;

   %fplot(@(m)the_eq2(m,area(10)),[1 5] ); grid on

   for i=1:length(area)
      chocked.mach(i)=fzero(@(m)the_eq(m,area(i)),[1 25]);
   endfor

% Pressure is the ratio p/P_stag
   chocked.pressure = ( 1 + 0.5*gamm1*(chocked.mach).^2 ).^(-gam/gamm1);
   %plot(xnozzle,pressure); grid on  
endfunction