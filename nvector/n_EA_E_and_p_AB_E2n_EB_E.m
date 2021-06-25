function [n_EB_E,z_EB] = n_EA_E_and_p_AB_E2n_EB_E(n_EA_E,p_AB_E,z_EA,a,f)
% n_EA_E_and_p_AB_E2n_EB_E From position A and delta, finds position B.
%    n_EB_E      = n_EA_E_and_p_AB_E2n_EB_E(n_EA_E,p_AB_E)
%   [n_EB_E,z_EB] = n_EA_E_and_p_AB_E2n_EB_E(n_EA_E,p_AB_E)
%   The n-vector for position A (n_EA_E) and the position-vector from position 
%   A to position B (p_AB_E) are given. The output is the n-vector of position 
%   B (n_EB_E) and depth of B (z_EB).
%   The calculation is exact, taking the ellipsity of the Earth into account.
%   It is also nonsingular as both n-vector and p-vector are nonsingular
%   (except for the center of the Earth).
%   The default ellipsoid model used is WGS-84, but other ellipsoids (or spheres)
%   might be specified.
%
%   [n_EB_E,z_EB] = n_EA_E_and_p_AB_E2n_EB_E(n_EA_E,p_AB_E,z_EA) Depth of A, z_EA, 
%   is also specified, z_EA = 0 is used when not spefified.
%
%   [n_EB_E,z_EB] = n_EA_E_and_p_AB_E2n_EB_E(n_EA_E,p_AB_E,z_EA,a)
%   Spherical Earth with radius a is used instead of WGS-84.
%
%   [n_EB_E,z_EB] = n_EA_E_and_p_AB_E2n_EB_E(n_EA_E,p_AB_E,z_EA,a,f)
%   Ellipsoidal Earth model with semi-major axis a and flattening f is used 
%   instead of WGS-84.
%
%   IN: 
%   n_EA_E:  [no unit] n-vector of position A, decomposed in E (3x1 vector).
%   p_AB_E:  [m]       Position vector from A to B, decomposed in E (3x1 vector).
%   z_EA:    [m]       (Optional, assumed to be zero if not given) Depth of system A, 
%                      relative to the ellipsoid (z_EA = -height).
%   a:       [m]       (Optional) Semi-major axis of the Earth ellipsoid
%   f:       [no unit] (Optional) Flattening of the Earth ellipsoid
%
%   OUT:
%   n_EB_E:  [no unit] n-vector of position B, decomposed in E (3x1 vector).
%   z_EB:    [m]       Depth of system B, relative to the ellipsoid (z_EB = -height).
%
%   The function also accepts vectorized form, i.e. n_EA_E and p_AB_E are 3xn matrixes, 
%   z_EA and z_EB are 1xn vectors and n_EB_E is a 3xn matrix.
%
%   See also n_EA_E_and_n_EB_E2p_AB_E, p_EB_E2n_EB_E, n_EB_E2p_EB_E.

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


% Originated: 2004.07.07 Kenneth Gade, FFI
% Modified:


%% INPUT HANDLING:

% The optional inputs a and f are forwarded to the kernel functions (which
% use the same syntax):
if nargin == 5 
    arg_a_f = {a,f};
elseif nargin == 4 
    arg_a_f = {a};
else
    arg_a_f = cell(0);
    if nargin == 2 % depth of A is not specified, using zero
        z_EA = zeros(1,size(n_EA_E,2));
    end
end

%% CALCULATIONS:

% Function 2. in Section 5.4 in Gade (2010):
p_EA_E = n_EB_E2p_EB_E(n_EA_E,z_EA,arg_a_f{:});
p_EB_E = p_EA_E+p_AB_E;
[n_EB_E,z_EB] = p_EB_E2n_EB_E(p_EB_E,arg_a_f{:});
