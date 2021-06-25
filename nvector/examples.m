function examples()
%   This file contains solutions to the 10 examples given at
%   www.navlab.net/nvector

%   The content of this file is based on the following publication:
%
%   Gade, K. (2010). A Nonsingular Horizontal Position Representation, The Journal
%   of Navigation, Volume 63, Issue 03, pp 395-417, July 2010.
%   (www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)

%   Copyright (c) 2021, Norwegian Defence Research Establishment (FFI)
%   All rights reserved.

%   Originated: 2015.03.26 Kenneth Gade, FFI

% NOTES:
% - All angles are by default assumed to be in radians, if an angle is
%   in degrees, the variable name has the following ending: _deg
%
% - The dot product (inner product) of vectors x and y is written dot(x,y)
%   here to make the code more readable for those unfamiliar with
%   Matlab. In Matlab one would normally write x'*y (i.e. x transposed
%   times y)


example1
example2
example3
example4
example5
example6
example7
example8
example9
example10

end


function example1() % Example 1: "A and B to delta"

% Positions A and B are given in (decimal) degrees and depths:

% Position A:
lat_EA_deg = 1;
long_EA_deg = 2;
z_EA = 3;

% Position B:
lat_EB_deg = 4;
long_EB_deg = 5;
z_EB = 6;

% Find the exact vector between the two positions, given in meters north,
% east, and down, i.e. find p_AB_N.


% SOLUTION:

% Step1: Convert to n-vectors (rad() converts to radians):
n_EA_E = lat_long2n_E(rad(lat_EA_deg),rad(long_EA_deg));
n_EB_E = lat_long2n_E(rad(lat_EB_deg),rad(long_EB_deg));

% Step2: Find p_AB_E (delta decomposed in E). WGS-84 ellipsoid is default:
[p_AB_E] = n_EA_E_and_n_EB_E2p_AB_E(n_EA_E,n_EB_E,z_EA,z_EB);

% Step3: Find R_EN for position A:
R_EN = n_E2R_EN(n_EA_E);

% Step4: Find p_AB_N
p_AB_N = R_EN'*p_AB_E;
% (Note the transpose of R_EN: The "closest-rule" says that when
% decomposing, the frame in the subscript of the rotation matrix that is
% closest to the vector, should equal the frame where the vector is
% decomposed. Thus the calculation R_NE*p_AB_E is correct, since the vector
% is decomposed in E, and E is closest to the vector. In the above example
% we only had R_EN, and thus we must transpose it: R_EN' = R_NE)

% Step5: Also find the direction (azimuth) to B, relative to north:
azimuth = atan2(p_AB_N(2),p_AB_N(1)); % positive angle about down-axis

disp(['Ex1: Delta north, east, down = ',num2str(p_AB_N(1)),', ', num2str(p_AB_N(2)),', ',num2str(p_AB_N(3)),' m']);
disp(['Ex1: Azimuth = ', num2str(deg(azimuth)),' deg']);
disp(' ');

end


function example2() % Example 2: "B and delta to C"

% delta vector from B to C, decomposed in B is given:
p_BC_B = [3000 2000 100]';

% Position and orientation of B is given:
n_EB_E = unit([1 2 3]'); % unit to get unit length of vector
z_EB = -400;
R_NB = zyx2R(rad(10),rad(20),rad(30)); % the three angles are yaw, pitch, and roll

% A custom reference ellipsoid is given (replacing WGS-84):
a = 6378135; f = 1/298.26; % (WGS-72)

% Find the position of C.


% SOLUTION:

% Step1: Find R_EN:
R_EN = n_E2R_EN(n_EB_E);

% Step2: Find R_EB, from R_EN and R_NB:
R_EB = R_EN*R_NB; % Note: closest frames cancel

% Step3: Decompose the delta vector in E:
p_BC_E = R_EB*p_BC_B; % no transpose of R_EB, since the vector is in B

% Step4: Find the position of C, using the functions that goes from one
% position and a delta, to a new position:
[n_EC_E,z_EC] = n_EA_E_and_p_AB_E2n_EB_E(n_EB_E,p_BC_E,z_EB,a,f);

% When displaying the resulting position for humans, it is more convenient
% to see lat, long:
[lat_EC,long_EC] = n_E2lat_long(n_EC_E);
% Here we also assume that the user wants the output to be height (= - depth):
disp(['Ex2: Pos C: lat, long = ',num2str(deg(lat_EC)),', ',num2str(deg(long_EC)),' deg, height = ',num2str(-z_EC),' m']);
disp(' ');

end


function example3() % Example 3: ECEF-vector to geodetic latitude

% Position B is given as p_EB_E ("ECEF-vector")

p_EB_E = 6371e3*[0.71 -0.72 0.1]'; % m

% Find position B as geodetic latitude, longitude and height

% SOLUTION:

% Find n-vector from the p-vector:
[n_EB_E,z_EB] = p_EB_E2n_EB_E(p_EB_E);

% Convert to lat, long and height:

[lat_EB,long_EB] = n_E2lat_long(n_EB_E);
h_EB = -z_EB;

disp(['Ex3: Pos B: lat, long = ',num2str(deg(lat_EB)),', ',num2str(deg(long_EB)),' deg, height = ',num2str(h_EB),' m']);
disp(' ');

end

function example4() % Example 4: geodetic latitude to ECEF-vector

% Position B is given with lat, long and height:
lat_EB_deg = 1;
long_EB_deg = 2;
h_EB = 3;

% Find the vector p_EB_E ("ECEF-vector")


% SOLUTION:

% Step1: Convert to n-vector:
n_EB_E = lat_long2n_E(rad(lat_EB_deg),rad(long_EB_deg));

% Step2: Find the ECEF-vector p_EB_E:
p_EB_E = n_EB_E2p_EB_E(n_EB_E,-h_EB);

disp(['Ex4: p_EB_E = [',num2str(p_EB_E(1)),', ',num2str(p_EB_E(2)),', ',num2str(p_EB_E(3)),']'' m']);
disp(' '); clear;

end


function example5() % Example 5: Surface distance

% Position A and B are given as n_EA_E and n_EB_E:
% Enter elements directly:
% n_EA_E = unit([1 0 -2]');
% n_EB_E = unit([-1 -2 0]');

% or input as lat/long in deg:
n_EA_E = lat_long2n_E(rad(88),rad(0));
n_EB_E = lat_long2n_E(rad(89),rad(-170));

r_Earth = 6371e3; % m, mean Earth radius

% SOLUTION:

% The great circle distance is given by equation (16) in Gade (2010):
% Well conditioned for all angles:
s_AB = (atan2(norm(cross(n_EA_E,n_EB_E)),dot(n_EA_E,n_EB_E))*r_Earth);

% % ill conditioned for small angels:
% s_AB_version1 = acos(dot(n_EA_E,n_EB_E))*r_Earth;
%
% % ill-conditioned for angles near pi/2 (and not valid above pi/2)
% s_AB_version2 = asin(norm(cross(n_EA_E,n_EB_E)))*r_Earth;

% The Euclidean distance is given by:
d_AB = norm(n_EB_E-n_EA_E)*r_Earth;

disp(['Ex5: Great circle distance = ',num2str(s_AB/1000),' km, Euclidean distance = ',num2str(d_AB/1000),' km']);
disp(' ');

end


function example6() % Example 6: Interpolated position

% Position B is given at time t0 as n_EB_E_t0 and at time t1 as n_EB_E_t1:
% Enter elements directly:
% n_EB_E_t0 = unit([1 0 -2]');
% n_EB_E_t1 = unit([-1 -2 0]');

% or input as lat/long in deg:
n_EB_E_t0 = lat_long2n_E(rad(89.9),rad(-150));
n_EB_E_t1 = lat_long2n_E(rad(89.9),rad(150));

% The times are given as:
t0 = 10;
t1 = 20;
ti = 16; % time of interpolation

% Find the interpolated position at time ti, n_EB_E_ti

% SOLUTION:

% Using standard interpolation:
n_EB_E_ti = unit(n_EB_E_t0+(ti-t0)*(n_EB_E_t1-n_EB_E_t0)/(t1-t0));

% When displaying the resulting position for humans, it is more convenient
% to see lat, long:
[lat_EB_ti,long_EB_ti] = n_E2lat_long(n_EB_E_ti);
disp(['Ex6: Interpolated position: lat, long = ',num2str(deg(lat_EB_ti)),', ',num2str(deg(long_EB_ti)),' deg']);
disp(' ');

end


function example7() % Example 7: Mean position/centre

% Three positions A, B and C are given:
% Enter elements directly:
% n_EA_E = unit([1 0 -2]');
% n_EB_E = unit([-1 -2 0]');
% n_EC_E = unit([0 -2 3]');

% or input as lat/long in degrees:
n_EA_E = lat_long2n_E(rad(90),rad(0));
n_EB_E = lat_long2n_E(rad(60),rad(10));
n_EC_E = lat_long2n_E(rad(50),rad(-20));

% Find the horizontal mean position, M:
n_EM_E = unit(n_EA_E+n_EB_E+n_EC_E);

% The result is best viewed with a figure that shows the n-vectors relative
% to an Earth-model:
disp('Ex7: See figure'); disp(' ');
plot_Earth_figure(n_EA_E,n_EB_E,n_EC_E,n_EM_E)

end


function example8() % Example 8: Position A and azimuth&distance to B

% Position A is given as n_EA_E:
% Enter elements directly:
% n_EA_E = unit([1 0 -2]');

% or input as lat/long in deg:
n_EA_E = lat_long2n_E(rad(80),rad(-90));

% The initial azimuth and great circle distance (s_AB), and Earth radius
% (r_Earth) are also given:
azimuth = rad(200);
s_AB = 1000; % m
r_Earth = 6371e3; % m, mean Earth radius

% Find the destination point B, as n_EB_E ("The direct/first geodetic
% problem" for a sphere)

% SOLUTION:

% Step1: Find unit vectors for north and east (see equations (9) and (10)
% in Gade (2010):
k_east_E = unit(cross(R_Ee'*[1 0 0]',n_EA_E));
k_north_E = cross(n_EA_E,k_east_E);

% Step2: Find the initial direction vector d_E:
d_E = k_north_E*cos(azimuth)+k_east_E*sin(azimuth);

% Step3: Find n_EB_E:
n_EB_E = n_EA_E*cos(s_AB/r_Earth)+d_E*sin(s_AB/r_Earth);

% When displaying the resulting position for humans, it is more convenient
% to see lat, long:
[lat_EB,long_EB] = n_E2lat_long(n_EB_E);
disp(['Ex8: Destination: lat, long = ',num2str(deg(lat_EB)),', ',num2str(deg(long_EB)),' deg']);
disp(' ');

end


function example9() % Example 9: Intersection

% Two paths A and B are given by two pairs of positions:
% Enter elements directly:
% n_EA1_E = unit([0 0 1]');
% n_EA2_E = unit([-1 0 1]');
% n_EB1_E = unit([-2 -2 4]');
% n_EB2_E = unit([-2 2 2]');


% or input as lat/long in deg:
n_EA1_E = lat_long2n_E(rad(50),rad(180));
n_EA2_E = lat_long2n_E(rad(90),rad(180));
n_EB1_E = lat_long2n_E(rad(60),rad(160));
n_EB2_E = lat_long2n_E(rad(80),rad(-140));

% Find the intersection between the two paths, n_EC_E:
n_EC_E_tmp = unit(cross(cross(n_EA1_E,n_EA2_E),cross(n_EB1_E,n_EB2_E)));

% n_EC_E_tmp is one of two solutions, the other is -n_EC_E_tmp. Select the
% one that is closest to n_EA1_E, by selecting sign from the dot product
% between n_EC_E_tmp and n_EA1_E:
n_EC_E = sign(dot(n_EC_E_tmp,n_EA1_E))*n_EC_E_tmp;

% When displaying the resulting position for humans, it is more convenient
% to see lat, long:
[lat_EC,long_EC] = n_E2lat_long(n_EC_E);
disp(['Ex9: Intersection: lat, long = ',num2str(deg(lat_EC)),', ',num2str(deg(long_EC)),' deg']);
disp(' ');

end


function example10() % Example 10: Cross track distance

% Position A1 and A2 and B are given as n_EA1_E, n_EA2_E, and n_EB_E:
% Enter elements directly:
% n_EA1_E = unit([1 0 -2]');
% n_EA2_E = unit([-1 -2 0]');
% n_EB_E = unit([0 -2 3]');

% or input as lat/long in deg:
n_EA1_E = lat_long2n_E(rad(0),rad(0));
n_EA2_E = lat_long2n_E(rad(10),rad(0));
n_EB_E = lat_long2n_E(rad(1),rad(0.1));

r_Earth = 6371e3; % m, mean Earth radius

% Find the cross track distance from path A to position B.

% SOLUTION:

% Find the unit normal to the great circle between n_EA1_E and n_EA2_E:
c_E = unit(cross(n_EA1_E,n_EA2_E));

% Find the great circle cross track distance: (acos(x) - pi/2 = -asin(x))
s_xt = -asin(dot(c_E,n_EB_E))*r_Earth;

% Find the Euclidean cross track distance:
d_xt = -dot(c_E,n_EB_E)*r_Earth;

disp(['Ex10: Cross track distance = ',num2str(s_xt),' m, Euclidean = ',num2str(d_xt),' m']);
disp(' ');

end


%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% The rest of the code is to plot the Earth with arrows for Example 7

function plot_Earth_figure(n_EA_E,n_EB_E,n_EC_E,n_EM_E)
% Plotting the Earth figure for Example 7.

Earth_radius_for_plotting = 0.7; % plot an Earth sphere with radius smaller
% than 1 (n-vector) to make the tip of the n-vectors visible.

% R_Ee selects correct E-axes, see R_Ee.m for details:
n_EA_E = R_Ee*n_EA_E; n_EB_E = R_Ee*n_EB_E; n_EC_E = R_Ee*n_EC_E; n_EM_E = R_Ee*n_EM_E;

figure(1)
clf;
hold on;

% To plot 3D arrows the function arrow3d.m, available from Matlab file
% exchange, written by Moshe Lindner is used (a copy of the function is at
% the end of this file):
arrow3d([0 -n_EA_E(3)],[0 n_EA_E(2)],[0 n_EA_E(1)],0.9,0.02,0.05,'r'); % A, n_EA_E
arrow3d([0 -n_EB_E(3)],[0 n_EB_E(2)],[0 n_EB_E(1)],0.9,0.02,0.05,'r'); % B, n_EB_E
arrow3d([0 -n_EC_E(3)],[0 n_EC_E(2)],[0 n_EC_E(1)],0.9,0.02,0.05,'r'); % C, n_EC_E
arrow3d([0 -n_EM_E(3)],[0 n_EM_E(2)],[0 n_EM_E(1)],0.9,0.02,0.05,'g'); % M, n_EM_E


%%%%%%%%%%%%%%%%%%%%%% Plotting a spherical Earth surface:
try
    % Loads a simple topographic model of Earth (available as part of default Matlab)
    load topo
    
    % Remove height information, only storing information about water or land:
    Earth_topo_binary = zeros(size(topo));
    Earth_topo_binary(topo > 0) = 1; % Set all positions above water to 1
    Earth_topo_binary(topo <= 0) = -1; % Set all positions below or equal to zero to -1
    Earth_topo_binary = [Earth_topo_binary(:,181:360) Earth_topo_binary(:,1:180)]; % Switch the halves to get correct mapping for our plot
    clear topo
catch
    Earth_topo_binary = zeros(2)-1;
end

% Number of elements in the sphere (in each direction)
n_of_Earth_plot_elements = 90;

% Data for a 3D Earth sphere
[Earth_plot_X,Earth_plot_Y,Earth_plot_Z] = sphere(n_of_Earth_plot_elements);
Earth_plot_X = Earth_plot_X*Earth_radius_for_plotting;
Earth_plot_Y = Earth_plot_Y*Earth_radius_for_plotting;
Earth_plot_Z = Earth_plot_Z*Earth_radius_for_plotting;

Earth_surface_properties.Cdata = Earth_topo_binary; % Color data
Earth_surface_properties.FaceLighting = 'gouraud'; % Smooth out light reflex
Earth_surface_properties.FaceColor = 'texture'; % Note that the 'facecolor' property needs to be set to
% 'texturemap' if the size of the z-data is different from the size of the data in the colormap (topo) that is loaded.
Earth_surface_properties.EdgeColor = 'none'; % Remove mesh
Earth_surface_properties.AmbientStrength = 0.1; % Ambient light strength
surface(Earth_plot_X,Earth_plot_Y,Earth_plot_Z,Earth_surface_properties);
material dull

colmap_Earth_binary = [0.6 0.6 1
    0.6 1 0.6]; % blue and green, coloring sea and land
colormap(colmap_Earth_binary)

% Add two light sources:
light('Position',[0 -1 0.5],'Style','infinite');
light('Position',[0 1 -0.5],'Style','infinite');

hold off;
grid on;
axis equal;
title('Example 7');
view(60,30);
rotate3d on;

end


% The following code is to plot 3D arrows. It is written by Moshe Lindner, and
% was retreieved from MATLAB Central File Exchange 2015.03.26
% https://www.mathworks.com/matlabcentral/fileexchange/28324-3d-arrow-plot

function [h] = arrow3d(x,y,z,head_frac,radii,radii2,colr)
%
% The function plotting 3-dimensional arrow
%
% h = arrow3d(x,y,z,head_frac,radii,radii2,colr)
%
% The inputs are:
%       x,y,z = vectors of the starting point and the ending point of the
%           arrow, e.g.:  x=[x_start, x_end]; y=[y_start, y_end];z=[z_start,z_end];
%       head_frac = fraction of the arrow length where the head should  start
%       radii = radius of the arrow
%       radii2 = radius of the arrow head (defult = radii*2)
%       colr =   color of the arrow, can be string of the color name, or RGB vector  (default='blue')
%
% The output is the handle of the surfaceplot graphics object.
% The settings of the plot can changed using: set(h, 'PropertyName', PropertyValue)
%
% example #1:
%        arrow3d([0 0],[0 0],[0 6],.5,3,4,[1 0 .5]);
% example #2:
%        arrow3d([2 0],[5 0],[0 -6],.2,3,5,'r');
% example #3:
%        h = arrow3d([1 0],[0 1],[-2 3],.8,3);
%        set(h,'facecolor',[1 0 0])
%
% Written by Moshe Lindner , Bar-Ilan University, Israel.
% July 2010 (C)

%
% Copyright (c) 2010, Moshe Lindner
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
%
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



if nargin==5
    radii2=radii*2;
    colr='blue';
elseif nargin==6
    colr='blue';
end
if size(x,1)==2
    x=x';
    y=y';
    z=z';
end

x(3)=x(2);
x(2)=x(1)+head_frac*(x(3)-x(1));
y(3)=y(2);
y(2)=y(1)+head_frac*(y(3)-y(1));
z(3)=z(2);
z(2)=z(1)+head_frac*(z(3)-z(1));
r=[x(1:2)',y(1:2)',z(1:2)'];

N=50;
dr=diff(r);
dr(end+1,:)=dr(end,:);
origin_shift=(ones(size(r))*(1+max(abs(r(:))))+[dr(:,1) 2*dr(:,2) -dr(:,3)]);
r=r+origin_shift;

normdr=(sqrt((dr(:,1).^2)+(dr(:,2).^2)+(dr(:,3).^2)));
normdr=[normdr,normdr,normdr];
dr=dr./normdr;
Pc=r;
n1=cross(dr,Pc);
normn1=(sqrt((n1(:,1).^2)+(n1(:,2).^2)+(n1(:,3).^2)));
normn1=[normn1,normn1,normn1];
n1=n1./normn1;
P1=n1+Pc;

X1=[];Y1=[];Z1=[];
j=1;
for theta=([0:N])*2*pi./(N)
    R1=Pc+radii*cos(theta).*(P1-Pc) + radii*sin(theta).*cross(dr,(P1-Pc)) -origin_shift;
    X1(2:3,j)=R1(:,1);
    Y1(2:3,j)=R1(:,2);
    Z1(2:3,j)=R1(:,3);
    j=j+1;
end

r=[x(2:3)',y(2:3)',z(2:3)'];

dr=diff(r);
dr(end+1,:)=dr(end,:);
origin_shift=(ones(size(r))*(1+max(abs(r(:))))+[dr(:,1) 2*dr(:,2) -dr(:,3)]);
r=r+origin_shift;

normdr=(sqrt((dr(:,1).^2)+(dr(:,2).^2)+(dr(:,3).^2)));
normdr=[normdr,normdr,normdr];
dr=dr./normdr;
Pc=r;
n1=cross(dr,Pc);
normn1=(sqrt((n1(:,1).^2)+(n1(:,2).^2)+(n1(:,3).^2)));
normn1=[normn1,normn1,normn1];
n1=n1./normn1;
P1=n1+Pc;

j=1;
for theta=([0:N])*2*pi./(N)
    R1=Pc+radii2*cos(theta).*(P1-Pc) + radii2*sin(theta).*cross(dr,(P1-Pc)) -origin_shift;
    X1(4:5,j)=R1(:,1);
    Y1(4:5,j)=R1(:,2);
    Z1(4:5,j)=R1(:,3);
    j=j+1;
end

X1(1,:)=X1(1,:)*0 + x(1);
Y1(1,:)=Y1(1,:)*0 + y(1);
Z1(1,:)=Z1(1,:)*0 + z(1);
X1(5,:)=X1(5,:)*0 + x(3);
Y1(5,:)=Y1(5,:)*0 + y(3);
Z1(5,:)=Z1(5,:)*0 + z(3);

h=surf(X1,Y1,Z1,'facecolor',colr,'edgecolor','none');
%camlight
lighting phong

end


