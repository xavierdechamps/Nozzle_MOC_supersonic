clear
clc
close all

%% Description
% PROGRAM to compute the supersonic isentropic flow inside a nozzle of known shape.
% The method of characteristics is used for this purpose.
%   params gives the parameters for the flow properties
%   geom gives the parameters for the geometry and the mesh.
%%

% Parameters for the flow field
params.gamma   = 1.2;    % Specific heat ratio = Cp / Cv
params.R       = 320;    % [J/kg-K] Gas constant
params.P       = 70.e5;  % [Pa] Stagnation pressure
params.PRatio  = 2;      % Ratio of static pressure at exit lip point [>1 to have a Prandtl-Meyer expansion]
                         % PRatio = Pstatic_lip / Pstatic_ambiant
params.T       = 3000;   % [K] Stagnation temperature

                       
% Parameters for the geometry
geom.delta = 1 ;         % [0/1] 0: planar nozzle
                         %       1: axisymmetric nozzle
geom.yt    = 1.  ;           % [m] Throat radius -> used as reference length for the whole nozzle
geom.rhou  = 2   * geom.yt ; % [m] Throat upstream radius of circular arc, required by MOC_2D_steady_irrotational_IVLINE
geom.rhod  = 0.5 * geom.yt ; % [m] Throat downstream radius of circular arc
geom.xe    = 10  * geom.yt ; % [m] Nozzle length
geom.xplume= 20   * geom.yt; % [m] Axial length of the plume from the nozzle exit lip point to the end (on the free pressure boundary)
geom.ta    = 15 ;            % [deg] Attachment angle between circular arc and line
geom.te    = 15 ;            % [deg] Exit lip point angle

% Parameters for the discretization
geom.NI            = 11 ;          % Number of points on the initial-value line
geom.circdownTheta = 1:1:geom.ta ; % Discretization of the circular arc downstream of the throat
geom.NIexpansion   = 10  ;          % Number of expansion waves, if any exists

% Parameters for the post-processing of the results
plots.patches      = 1; % [0/1] Plot the 2D patches to visualize the characteristics and 
                        %       the solution inside the nozzle
plots.patches_data = 1; % 0: plot only the outlines of the patches -> no colour
                        % 1: plot the Mach number 
                        % 2:      the static pressure
                        % 3:      the static temperature
                        % 4:      the density
                        % 5:      the amplitude of the velocity
                        % 6:      flow direction [deg] in comparison with the axial direction
plots.patches_xlim = 40 * geom.yt; % Abscissa above which the patch plot is cut

addpath('./src/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry of the downstream circular arc
geom.circdownX =           geom.rhod*   sind( geom.circdownTheta ) ;
geom.circdownY = geom.yt + geom.rhod*(1-cosd( geom.circdownTheta)) ;

% The initial-value line is chosen as the line where V=0
ythroat = 0:geom.yt/(geom.NI-1):geom.yt;
% xsonic contains the absissae for the sonic line Mach=1 at the throat, not used here
% xvnull contains the absissae the the line V=0 at the throat
% uvnull contains the value of U on the line V=0
[xsonic,xvnull,uvnull] = MOC_2D_steady_irrotational_IVLINE ( geom , params , ythroat );
vvnull = 0;

ythroat(1) += 1.e-6; % To avoid singularity on axis
uvnull(1)  += 1.e-6; % To avoid singularity on axis

X(1,1) = xvnull(1);
Y(1,1) = ythroat(1);
U(1,1) = uvnull(1);
V(1,1) = vvnull;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  E X T E N D    T H E    I N I T I A L    V A L U E - L I N E  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extend the solution from the initial-value line, fill-in the matrices X, Y, U and V
% which contains the (x,y) coordinates of the intersections of the characteristics
% and the (u,v) velocity components at these points
%
disp('Computing the extent from the initial-value line...')
[indI,X,Y,U,V,LENG_INDI,ind_patch_tri,ind_patch_quad,patch_i_tri,patch_j_tri,patch_i_quad,patch_j_quad] = ...
                                    MOC_2D_steady_irrotational_extend_initial_line(geom,params,...
                                                                                   xvnull,ythroat,uvnull,vvnull,
                                                                                   X,Y,U,V) ;
LOC_PATCHES(1,:) = [ind_patch_tri  ind_patch_quad]; % Used to locate the different regions when plotting the patches

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  E X T E N D    T H E    I N I T I A L    E X P A N S I O N  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse wall point method for the downwards circular arc with prespecified points
% which contains the (x,y) coordinates of the intersections of the characteristics
% and the (u,v) velocity components at these points
%
disp('Computing the extent of the flow field determined by the initial-expansion contour...')
[indI,X,Y,U,V,LENG_INDI,ind_patch_tri,ind_patch_quad,patch_i_tri,patch_j_tri,patch_i_quad,patch_j_quad] = ...
                                    MOC_2D_steady_irrotational_extend_initial_expansion(geom,params,indI,LENG_INDI,...
                                                                                        X,Y,U,V,...
                                                                                        ind_patch_tri,ind_patch_quad,...
                                                                                        patch_i_tri,patch_j_tri,patch_i_quad,patch_j_quad);
LOC_PATCHES(2,:) = [ind_patch_tri  ind_patch_quad];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  E X T E N D    T H E    F I N A L    E X P A N S I O N  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct wall point method for the linear part of the nozzle downstream of the circular arc
% which contains the (x,y) coordinates of the intersections of the characteristics
% and the (u,v) velocity components at these points
%
disp('Computing the extent of the flow field determined by the wall...')
[indI,X,Y,U,V,LENG_INDI,ind_patch_tri,ind_patch_quad,patch_i_tri,patch_j_tri,patch_i_quad,patch_j_quad] = ...
                                    MOC_2D_steady_irrotational_extend_expansion(geom,params,indI,LENG_INDI,...
                                                                                X,Y,U,V,...
                                                                                ind_patch_tri,ind_patch_quad,...
                                                                                patch_i_tri,patch_j_tri,patch_i_quad,patch_j_quad);
LOC_PATCHES(3,:) = [ind_patch_tri  ind_patch_quad];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  E X T E N D    T H E    P L U M E  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the static pressure at the nozzle lip point. If static pressure > ambiant static pressure: Prandtl-Meyer expansion
%                                                    If static pressure < ambiant static pressure: weak shock
% The Prandtl-Meyer expansion is discretized by geom.NIexpansion lines.
% The free-pressure point method is used to obtain the shape of the plume.
%
[indI,X,Y,U,V,LENG_INDI,ind_patch_tri,ind_patch_quad,patch_i_tri,patch_j_tri,patch_i_quad,patch_j_quad] = ...
                                    MOC_2D_steady_irrotational_extend_plume(geom,params,indI,LENG_INDI,...
                                                                            X,Y,U,V,...
                                                                            ind_patch_tri,ind_patch_quad,...
                                                                            patch_i_tri,patch_j_tri,patch_i_quad,patch_j_quad);
LOC_PATCHES(4,:) = [ind_patch_tri  ind_patch_quad];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  P O S T - P R O C E S S    T H E    R E S U L T S  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Plotting the results...')

%[Mach,pressure,density,temperature,sound] = MOC_2D_steady_irrotational_get_thermo(U,V,params);

[fig1,fig2,fig3] = MOC_2D_steady_irrotational_postprocess(geom,params,plots,indI,LENG_INDI,...
                                                          X,Y,U,V,...
                                                          ind_patch_tri,ind_patch_quad,...
                                                          patch_i_tri,patch_j_tri,patch_i_quad,patch_j_quad,LOC_PATCHES);

%saveas(fig1,'Characteristics.pdf') % pdf / jpg
%saveas(fig2,'Pressure.pdf')        % pdf / jpg
%saveas(fig3,'Mach_number.pdf')     % pdf / jpg