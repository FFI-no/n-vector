function [n_EB_E,z_EB] = p_EB_E2n_EB_E(p_EB_E,a,f)
% p_EB_E2n_EB_E Converts Cartesian position vector in meters to n-vector.
%   [n_EB_E,z_EB] = p_EB_E2n_EB_E(p_EB_E)
%   The position of B (typically body) relative to E (typically Earth) is
%   given into this function as cartesian position vector p_EB_E, in meters
%   ("ECEF-vector"). The function converts to n-vector, n_EB_E and its
%   depth, z_EB.
%   The calculation is exact, taking the ellipsity of the Earth into account.
%   It is also nonsingular as both n-vector and p-vector are nonsingular
%   (except for the center of the Earth).
%   The default ellipsoid model used is WGS-84, but other ellipsoids (or spheres)
%   might be specified.
%
%   [n_EB_E,z_EB] = p_EB_E2n_EB_E(p_EB_E,a) Spherical Earth with radius a is
%   used instead of WGS-84.
%
%   [n_EB_E,z_EB] = p_EB_E2n_EB_E(p_EB_E,a,f) Ellipsoidal Earth model with
%   semi-major axis a and flattening f is used instead of WGS-84.
%
%
%   IN:
%   p_EB_E: [m]       Cartesian position vector from E to B, decomposed in E (3x1 vector).
%   a:      [m]       (Optional) Semi-major axis of the Earth ellipsoid
%   f:      [no unit] (Optional) Flattening of the Earth ellipsoid
%
%   OUT:
%   n_EB_E: [no unit] n-vector  representation of position B, decomposed in E (3x1 vector).
%   z_EB:   [m]       Depth of system B relative to the ellipsoid (z_EB = -height).
%
%
%   The function also accepts vectorized form, i.e. p_EB_E is a 3xn matrix,
%   n_EB_E is a 3xn matrix and z_EB is a 1xn vector.
%
%   See also n_EB_E2p_EB_E, n_EA_E_and_p_AB_E2n_EB_E, n_EA_E_and_n_EB_E2p_AB_E.

% MIT License:
%
% Copyright (c) 2021, Norwegian Defence Research Establishment (FFI)
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.


%   This file is part of NavLab and is available from www.navlab.net/nvector
%
%   The content of this file is based on the following publication:
%
%   Gade, K. (2010). A Nonsingular Horizontal Position Representation, The Journal 
%   of Navigation, Volume 63, Issue 03, pp 395-417, July 2010. 
%   (www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)
%
%   This paper should be cited in publications using this file.


% Originated: 2004.11.17 Kenneth Gade and Brita Hafskjold, FFI
% Modified:   2007.03.02 Brita Hafskjold Gade, FFI
%                Replaced formulas to get full numerical accuracy at all positions
% Modified:   2014.08.22 Kenneth Gade, FFI:
%                Added possibility of vectorized input/output


%% INPUT HANDLING:

p_EB_E = R_Ee*p_EB_E; % R_Ee selects correct E-axes, see R_Ee.m for details

if nargin == 1 % WGS-84 ellipsoid is used
    a = 6378137; % the equatorial radius of the Earth-ellipsoid
    f = 1/298.257223563; % the flattening of the Earth-ellipsoid
    
elseif nargin == 2% custom sphere is specified:
    f = 0;
else % custom ellipsoid is specified:
    % Backward compatibility:
    % Previously, custom ellipsoid was specified by a and b in this function.
    % However, for more spherical globes than the Earth, or if f has more
    % decimals than in WGS-84, using f and a as input will give better
    % numerical precision than a and b.
    
    % old input number 3, 4: Polar_semi_axis (b), equatorial_semi_axis (a)
    if f > 1e6 % Checks if a is given as f (=old input)
        f_new = 1-a/f; a = f; f = f_new; % switch old inputs to new format
    end
end


%% CALCULATIONS:

%e_2 = eccentricity^2
e_2 = 2*f-f^2; % = 1-b^2/a^2;


% The following code implements equation (23) from Gade (2010):

R_2 = p_EB_E(2,:).^2 + p_EB_E(3,:).^2;
R = sqrt(R_2); % R = component of p_EB_E in the equatorial plane

p = R_2/a^2;
q = (1-e_2)/(a^2) * p_EB_E(1,:).^2;
r = (p+q-e_2^2)/6;

s = e_2^2*p.*q./(4*r.^3);
t = nthroot((1 + s + sqrt(s.*(2+s))),3);
u = r.*(1 + t + 1./t);
v = sqrt(u.^2 + e_2^2 .* q);

w = e_2*(u + v - q)./(2*v);
k = sqrt(u + v + w.^2) - w;
d = k.*R./(k + e_2);

% Calculate height:
height = (k + e_2 - 1)./k  .* sqrt(d.^2 + p_EB_E(1,:).^2);

temp = 1./sqrt(d.^2 + p_EB_E(1,:).^2);

n_EB_E_x = temp .* p_EB_E(1,:);
n_EB_E_y = temp .* k./(k+e_2) .* p_EB_E(2,:);
n_EB_E_z = temp .* k./(k+e_2) .* p_EB_E(3,:);

n_EB_E = [n_EB_E_x;n_EB_E_y;n_EB_E_z];

% Ensure unit length:
n_EB_E = unit(R_Ee'*n_EB_E);

z_EB = -height;
