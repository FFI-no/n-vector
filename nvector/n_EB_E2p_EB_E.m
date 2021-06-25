function p_EB_E = n_EB_E2p_EB_E(n_EB_E,z_EB,a,f)
% n_EB_E2p_EB_E Converts n-vector to Cartesian position vector in meters.
%   p_EB_E = n_EB_E2p_EB_E(n_EB_E)
%   The position of B (typically body) relative to E (typically Earth) is
%   given into this function as n-vector, n_EB_E. The function converts
%   to cartesian position vector ("ECEF-vector"), p_EB_E, in meters.
%   The calculation is exact, taking the ellipsity of the Earth into account.
%   It is also nonsingular as both n-vector and p-vector are nonsingular
%   (except for the center of the Earth).
%   The default ellipsoid model used is WGS-84, but other ellipsoids (or spheres)
%   might be specified.
%
%   p_EB_E = n_EB_E2p_EB_E(n_EB_E,z_EB) Depth of B, z_EB, is also specified,
%   z_EB = 0 is used when not specified.
%
%   p_EB_E = n_EB_E2p_EB_E(n_EB_E,z_EB,a) Spherical Earth with radius a is
%   used instead of WGS-84.
%
%   p_EB_E = n_EB_E2p_EB_E(n_EB_E,z_EB,a,f) Ellipsoidal Earth model with
%   semi-major axis a and flattening f is used instead of WGS-84.
%
%
%   IN:
%   n_EB_E:  [no unit] n-vector of position B, decomposed in E (3x1 vector).
%   z_EB:    [m]       (Optional, assumed to be zero if not given) Depth of system B,
%                      relative to the ellipsoid (z_EB = -height)
%   a:       [m]       (Optional) Semi-major axis of the Earth ellipsoid
%   f:       [no unit] (Optional) Flattening of the Earth ellipsoid
%
%   OUT:
%   p_EB_E:  [m]       Cartesian position vector from E to B, decomposed in E (3x1 vector).
%
%   The function also accepts vectorized form, i.e. n_EB_E is a 3xn matrix, z_EB is
%   a 1xn vector and p_EB_E is a 3xn matrix.
%
%   See also p_EB_E2n_EB_E, n_EA_E_and_p_AB_E2n_EB_E, n_EA_E_and_n_EB_E2p_AB_E.

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
% Modified:   2015.02.20 Kenneth Gade, FFI: Added possibility of using alternative E axes


%% INPUT HANDLING:

length_deviation_warning_limit = 0.1; % n-vector should have length = 1
% i.e. norm(n_E) = 1. A deviation from 1 exceeding this limit gives a warning.
% This function uses unit() to ensure length = 1, thus the warning is included
% only to give a notice in cases where a wrong input is given unintentionally 
% (i.e. the input is not even approximately a unit vector).

length_deviation = abs(norm(n_EB_E(:,1))-1); % If a matrix of n-vectors is input, 
% only first is controlled to save time (assuming advanced users input correct n-vectors)

if length_deviation>length_deviation_warning_limit
    warning('n_EB_E2p_EB_E: norm(n_EB_E)~= 1 ! Error is: %g',length_deviation); 
end

n_EB_E = unit(R_Ee*n_EB_E); % Ensures unit length. R_Ee selects correct E-axes, see R_Ee.m for details.
% Note: In code where the norm of the input n_EB_E is guaranteed to be 1,
% the use of the unit-function can be removed, to gain some speed.

if nargin <= 2 % WGS-84 ellipsoid is used
    if nargin == 1 % depth is not specified, using zero
        z_EB = zeros(1,size(n_EB_E,2));
    end
    a = 6378137; % the equatorial radius of the Earth-ellipsoid
    f = 1/298.257223563; % the flattening of the Earth-ellipsoid
    
    
elseif nargin == 3 %custom sphere is specified:  
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

% semi-minor axis:
b = a*(1-f);

% The following code implements equation (22) in Gade (2010):

denominator = sqrt(n_EB_E(1,:).^2 ...
    + n_EB_E(2,:).^2 / (1-f)^2 ...
    + n_EB_E(3,:).^2 / (1-f)^2);

% We first calculate the position at the origin of coordinate system L,
% which has the same n-vector as B (n_EL_E = n_EB_E),
% but lies at the surface of the Earth (z_EL = 0).

p_EL_E = [b./denominator.*n_EB_E(1,:);
    b./denominator.*n_EB_E(2,:)/(1-f)^2;
    b./denominator.*n_EB_E(3,:)/(1-f)^2];
% (The factor ellipsoid_semiaxis_Ex./denominator.* is put inside to make it work on vectorized form)

p_EB_E = R_Ee'*(p_EL_E-[n_EB_E(1,:).*z_EB
    n_EB_E(2,:).*z_EB
    n_EB_E(3,:).*z_EB]);
% If not vectorized form needed, use p_EB_E = p_EL_E - n_EB_E*z_EB;

