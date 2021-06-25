function n_E = lat_long2n_E(latitude,longitude)
% lat_long2n_E Converts latitude and longitude to n-vector.
%   n_E = lat_long2n_E(latitude,longitude) 
%   n-vector (n_E) is calculated from (geodetic) latitude and longitude.
%
%   IN:
%   latitude:  [rad]     Geodetic latitude
%   longitude: [rad]
%
%   OUT:
%   n_E:       [no unit] n-vector decomposed in E (3x1 vector)
%
%   The function also accepts vectors (1xm) with lat and long, then a 3xm 
%   matrix of n-vectors is returned.
%
%   See also n_E2lat_long.

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
%   Gade, K. (2010). A Nonsingular Horizontal Position Representation, The
%   Journal of Navigation, Volume 63, Issue 03, pp 395-417, July 2010.
%   (www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)
%
%   This paper should be cited in publications using this file.


% Originated: 1999.02.23 Kenneth Gade, FFI
% Modified:   2015.02.20 Kenneth Gade, FFI: Added possibility of using alternative E axes


% Equation (3) from Gade (2010):
n_E = R_Ee'*[sin(latitude) % R_Ee selects correct E-axes, see R_Ee.m for details
    sin(longitude).*cos(latitude)
    -cos(longitude).*cos(latitude)];
