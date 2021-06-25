function [latitude,longitude] = n_E2lat_long(n_E)
% n_E2lat_long Converts n-vector to latitude and longitude.
%   [latitude,longitude] = n_E2lat_long(n_E) 
%   Geodetic latitude and longitude are calculated from n-vector (n_E).
%
%   IN:
%   n_E:       [no unit] n-vector decomposed in E (3x1 vector)
%
%   OUT:
%   latitude:  [rad]     Geodetic latitude
%   longitude: [rad]
%
%   The function also accepts vectorized form (i.e. a 3xm matrix of n-vectors 
%   is input, returning 1xm vectors of latitude and longitude)
%
%   See also lat_long2n_E.

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
% Modified:   2004.11.23 Kenneth Gade, FFI: Accepts vectorized input
% Modified:   2015.02.20 Kenneth Gade, FFI: Added possibility of using alternative E axes


%% INPUT HANDLING:

length_deviation_warning_limit = 0.1; % n-vector should have length = 1
% i.e. norm(n_E) = 1. A deviation from 1 exceeding this limit gives a warning.
% This function only depends of the direction of n-vector, thus the warning 
% is included only to give a notice in cases where a wrong input is given 
% unintentionally (i.e. the input is not even approximately a unit vector).

length_deviation = abs(norm(n_E(:,1))-1); % If a matrix of n-vectors is input, 
% only the first n-vector is controlled to save time (assuming advanced
% users input correct n-vectors)

if length_deviation > length_deviation_warning_limit
    warning('n_E2lat_long: norm(n_E) ~= 1 ! Error is: %g',length_deviation); 
end

n_E = R_Ee*n_E; % R_Ee selects correct E-axes, see R_Ee.m for details


%% CALCULATIONS:

% Equation (5) in Gade (2010):
longitude = atan2(n_E(2,:),-n_E(3,:));

% Equation (6) in Gade (2010) (Robust numerical solution)
equatorial_component = sqrt(n_E(2,:).^2+n_E(3,:).^2); % vector component in the equatorial plane
latitude = atan2(n_E(1,:),equatorial_component); %atan() could also be used since latitude is within [-pi/2,pi/2]

% latitude = asin(n_E(1)) is a theoretical solution, but close to the Poles it
% is ill-conditioned which may lead to numerical inaccuracies (and it will give imaginary results for norm(n_E)>1)


