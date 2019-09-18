clear
clc
close all

format long e

addpath('.');
rad2deg = 180/pi;
gamma = 1.4;
Mach_out = 2.4 ;
num_rays = 7 ;
throat_rad = 1;

[nu_out,mu_out] = get_prandtl_meyer_function(Mach_out,gamma) ;
theta_max = 0.5 * nu_out ;

theta_init = theta_max - floor(theta_max);
theta(1,:) = theta_init : (theta_max-theta_init)/(num_rays-1) : theta_max ;
nu = theta ;
Mach = get_Mach_from_nu(nu,gamma);
[nu_dummy,mu] = get_prandtl_meyer_function(Mach,gamma) ;

Kminus = theta + nu;
Kplus  = theta - nu;

abs(1,1) = throat_rad / tand( abs(theta(1,1)-mu(1,1)) );
ord(1,1) = 0.;
for k=2:num_rays
  meanSlopeplus  = 0.5 * ( theta(1,k-1) + theta(1,k) + mu(1,k-1) + mu(1,k) ) / rad2deg ;
%  meanSlopeminus = 0.5 * ( theta(1,k-1) + theta(1,k) - mu(1,k-1) - mu(1,k) ) / rad2deg;
  meanSlopeminus = ( theta(1,k) - mu(1,k) ) / rad2deg ;
  
  abs(1,k) = (throat_rad-ord(1,k-1)+abs(1,k-1)*meanSlopeplus) / (meanSlopeplus-meanSlopeminus);
  ord(1,k) = throat_rad + meanSlopeminus*abs(1,k);
endfor

k = num_rays;
Wall(1,:) = [ Kminus(k) Kplus(k) theta(k) nu(k) Mach(k) mu(k) ] ;
meanSlopeplus = 0.5 * ( theta(1,k) + Wall(1,3) + mu(1,k) + Wall(1,6) ) / rad2deg ;
meanSlopeWall = 0.5 * ( theta_max + Wall(1,3) ) / rad2deg ;
absWall(1) = (throat_rad-ord(1,k)+abs(1,k)*meanSlopeplus) / (meanSlopeplus-meanSlopeWall);
ordWall(1) = throat_rad + meanSlopeWall*absWall(1);

for j=2:num_rays
  lk = (num_rays-j+1) ;
  for k=1:lk
    
    if (k==1)
      theta(j,k)         = 0. ;
      Kminus(j,k)        = Kminus(j-1,k+1) ;
      nu(j,k)            = Kminus(j,k) - theta(j,k) ;
      Kplus(j,k)         = theta(j,k)  - nu(j,k) ;
      Mach(j,k)          = get_Mach_from_nu(nu(j,k),gamma) ;
      [nu_dummy,mu(j,k)] = get_prandtl_meyer_function(Mach(j,k),gamma) ;
      
      meanSlopeminus     = 0.5 * ( theta(j-1,k+1) + theta(j,k) - mu(j-1,k+1) - mu(j,k) ) / rad2deg;
      abs(j,k)           = abs(j-1,k+1) - ord(j-1,k+1)/meanSlopeminus ;
      ord(j,k)           = 0. ;
    else
      Kplus(j,k)         = Kplus(j,k-1) ;
      Kminus(j,k)        = Kminus(j-1,k+1) ;
      theta(j,k)         = 0.5 * ( Kminus(j,k) + Kplus(j,k) );
      nu(j,k)            = 0.5 * ( Kminus(j,k) - Kplus(j,k) );
      Mach(j,k)          = get_Mach_from_nu(nu(j,k),gamma) ;
      [nu_dummy,mu(j,k)] = get_prandtl_meyer_function(Mach(j,k),gamma) ;
      
      meanSlopeplus      = 0.5 * ( theta(j,k-1)   + theta(j,k) + mu(j,k-1)   + mu(j,k) ) / rad2deg ;
      meanSlopeminus     = 0.5 * ( theta(j-1,k+1) + theta(j,k) - mu(j-1,k+1) - mu(j,k) ) / rad2deg ;
      abs(j,k)           = (ord(j-1,k+1)-ord(j,k-1)+abs(j,k-1)*meanSlopeplus-abs(j-1,k+1)*meanSlopeminus) / (meanSlopeplus-meanSlopeminus);
      ord(j,k)           = ord(j-1,k+1) + meanSlopeminus*(abs(j,k)-abs(j-1,k+1));
    endif
    
  endfor
  
  Wall(j,:) = [ Kminus(j,lk) Kplus(j,lk) theta(j,lk) nu(j,lk) Mach(j,lk) mu(j,lk) ] ;
  
  meanSlopeplus = 0.5 * ( theta(j,lk) + Wall(j,3) + mu(j,lk) + Wall(j,6) ) / rad2deg ;
  meanSlopeWall = 0.5 * ( Wall(j-1,3) + Wall(j,3) ) / rad2deg ;
  absWall(j) = (ordWall(j-1)-ord(j,lk)+abs(j,lk)*meanSlopeplus-absWall(j-1)*meanSlopeWall) / (meanSlopeplus-meanSlopeWall);
  ordWall(j) = ordWall(j-1) + meanSlopeWall*(absWall(j)-absWall(j-1));
endfor

ind=1; len = num_rays;
for j=1:num_rays
  absp(ind:ind+len-1) = abs(j,1:len);
  ordp(ind:ind+len-1) = ord(j,1:len);
  Machp(ind:ind+len-1) = Mach(j,1:len);
  ind=ind+len;
  len--;
endfor
%plot(abs,ord,'k+',[0 absWall],[throat_rad ordWall],'r+'); grid on; axis equal

scatter(absp,ordp,25,Machp,'filled'); hold on; grid on; colorbar; 
ylim([0 max(ordWall)]); 
plot([0 absWall],[throat_rad ordWall],'k')

area_ratio = (ordWall(end)/throat_rad)