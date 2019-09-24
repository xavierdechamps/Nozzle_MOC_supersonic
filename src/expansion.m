clear;
clc

M1 = 3.518326326701526;
theta1 = 15 ;
gamma = 1.2;

Pstag   = 150.e5; % [Pa] Stagnation pressure
Pstatic = 101325;

p1static = Pstag / ( 1. + 0.5*(gamma-1)*M1.^2 )^(gamma/(gamma-1)) ;

[nu1,mu1] = get_prandtl_meyer_function(M1,gamma);

% Initial values, to launch the iterative process
delta = 20; % deflection angle
% nu2 = nu1 + delta
[M2] = get_Mach_from_nu(nu1+delta,gamma);
p2static = Pstag / ( 1. + 0.5*(gamma-1)*M2.^2 )^(gamma/(gamma-1)) ;


FinalWaveExpansion(delta,gamma,nu1,Pstag,Pstatic) ;

delta = fzero(@(delta12) FinalWaveExpansion(delta12,gamma,nu1,Pstag,Pstatic),[0 delta],optimset('TolX',1e-12))