function R_EL = n_E_and_wa2R_EL(n_E,wander_azimuth)
% n_E_and_wa2R_EL Finds R_EL from n-vector and wander azimuth angle.
%   R_EL = n_E_and_wa2R_EL(n_E,wander_azimuth) 
%   Calculates the rotation matrix (direction cosine matrix) R_EL using
%   n-vector (n_E) and the wander azimuth angle.
%   When wander_azimuth = 0, we have that N = L (See Table 2 in Gade (2010) for
%   details)
%
%   IN: 
%   n_E:        [no unit] n-vector decomposed in E (3x1 vector)
%   wander_azimuth: [rad] The angle between L's x-axis and north, pos about L's z-axis
%
%   OUT:
%   R_EL:       [no unit] The resulting rotation matrix (3x3)
%
%   See also R_EL2n_E, R_EN2n_E, n_E2R_EN.

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


% Originated: 1999.02.23 Kenneth Gade, FFI
% Modified:   2015.02.20 Kenneth Gade, FFI: Added possibility of using alternative E axes


[latitude,longitude] = n_E2lat_long(n_E);

% Longitude, -latitude, and wander azimuth are the x-y-z Euler angles (about
% new axes) for R_EL. See also the second paragraph of Section 5.2 in Gade (2010):
R_EL = R_Ee'*xyz2R(longitude,-latitude,wander_azimuth); % R_Ee selects correct E-axes, see R_Ee.m for details


