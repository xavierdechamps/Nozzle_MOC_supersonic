%function oblique_shock()
 clear
 clc
 close all

 Mach = [2 3 5];
 gamma = 1.4;

 theta_min = 1*pi/180;
 theta_steps = 100;

 % Equations to solve
 equation_max_theta       = @(gamma,Mach,beta)       tan( beta ) .* ( 0.5 *( (gamma+cos(2*beta))*Mach^2 + 2) ) ./ ( (Mach*sin(beta)).^2 - 1 ) ;
 equation_Theta_Beta_Mach = @(gamma,Mach,theta,beta) 2 * cot( beta ) .* ( (Mach*sin(beta)).^2 - 1 ) ./ ( (gamma+cos(2*beta))*Mach^2 + 2) - tan(theta);
 
 Mach_min = 1.1;
 Mach_max = Mach(end);
 Mach_theta_max = Mach_min : (Mach_max-Mach_min)/100 : Mach_max ;
 k=1;
 for mach_j=Mach_theta_max
   beta_max = fminbnd ( @(beta) equation_max_theta(gamma,mach_j,beta) , pi/3 , pi/2-0.0001 ) ;
   theta_max = acot(equation_max_theta(gamma,mach_j,beta_max));
   machmax(k,:) = [ mach_j beta_max*180/pi theta_max*180/pi ] ;
   k++;
 end
 
 figure()
 j=1;
 for mach_j = Mach
  % Get theta_max and beta (rad) corresponding to theta_max
   beta_max = fminbnd ( @(beta) equation_max_theta(gamma,mach_j,beta) , pi/4 , pi/2-0.001 ) ;
   theta_max = acot(equation_max_theta(gamma,mach_j,beta_max));
   
  % Get the value of theta_max (rad)
   theta_max_search = theta_max - 0.00001;
   theta = theta_min : (theta_max_search-theta_min)/theta_steps : theta_max_search ;
   i=1;
   for theta_i = theta
    try
  %    Initial guess of Mach number
      guessBeta = [pi/180 beta_max];
      
      beta_r = guessBeta(1) : pi/100 : guessBeta(2);
      fct = equation_Theta_Beta_Mach(gamma,mach_j,theta_i,beta_r);
      plot(beta_r,fct); grid on
      
      betaEval_1 = fzero(@(beta) equation_Theta_Beta_Mach(gamma,mach_j,theta_i,beta), guessBeta);
      
      guessBeta = [betaEval_1+0.0001 pi/2];
      betaEval_2 = fzero(@(beta) equation_Theta_Beta_Mach(gamma,mach_j,theta_i,beta), guessBeta);
      
      betadev(j,i,:) = [ theta_i betaEval_1 betaEval_2 ]*180/pi ;
      i++;
    catch e
      disp('fzero failed')
      break;
    end
   end
   plot(betadev(j,1:(i-1),1),betadev(j,1:(i-1),2),'r',betadev(j,1:(i-1),1),betadev(j,1:(i-1),3),'b'); grid on;hold on; drawnow;
   j++;
 end
 
 plot(machmax(:,3),machmax(:,2),'k');
 title(['Mach number = ',num2str(Mach),', max theta = ',num2str(theta_max*180/pi)])
 xlabel('Shock wave angle beta (deg)') ; ylabel('Deflection angle theta (deg)') ; 