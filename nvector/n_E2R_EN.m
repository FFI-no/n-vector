function R_EN = n_E2R_EN(n_E)
% n_E2R_EN Finds the rotation matrix R_EN from n-vector.
%   R_EN = n_E2R_EN(n_E) 
%   The rotation matrix (direction cosine matrix) R_EN is calculated based
%   on n-vector (n_E).
%
%   IN:
%   n_E:   [no unit] n-vector decomposed in E (3x1 vector)
%
%   OUT:
%   R_EN:  [no unit] The resulting rotation matrix (direction cosine matrix)
%
%   See also R_EN2n_E, n_E_and_wa2R_EL, R_EL2n_E.

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


% Originated: 2015.02.23 Kenneth Gade, FFI
% Modified:


%% INPUT HANDLING:

length_deviation_warning_limit = 0.1; % n-vector should have length = 1
% i.e. norm(n_E) = 1. A deviation from 1 exceeding this limit gives a warning.
% This function uses unit() to ensure length = 1, thus the warning is included
% only to give a notice in cases where a wrong input is given unintentionally 
% (i.e. the input is not even approximately a unit vector)

length_deviation = abs(norm(n_E(:,1))-1); % If a matrix of n-vectors is input, 
% only the first n-vector is controlled to save time (assuming advanced
% users input correct n-vectors)

if length_deviation>length_deviation_warning_limit
    warning('n_E2R_EN: norm(n_E)~= 1 ! Error is: %g',length_deviation); 
end

n_E = unit(R_Ee*n_E); % Ensures unit length. R_Ee selects correct E-axes, see R_Ee.m for details.
% Note: In code where the norm of the input n_EB_E is guaranteed to be 1,
% the use of the unit-function can be removed, to gain some speed.


%% CALCULATIONS:

% N coordinate frame (North-East-Down) is defined in Table 2 in Gade (2010)

% R_EN is constructed by the following three column vectors: The x, y and z
% basis vectors (axes) of N, each decomposed in E.

% Find z-axis of N (Nz):
Nz_E = -n_E; % z-axis of N (down) points opposite to n-vector

% Find y-axis of N (East)(remember that N is singular at Poles)
% Equation (9) in Gade (2010):
Ny_E_direction = cross([1 0 0]',n_E); % Ny points perpendicular to the plane
% formed by n-vector and Earth's spin axis                                   
if norm(Ny_E_direction) ~= 0 % outside Poles:    
    Ny_E = unit(Ny_E_direction);
else % Pole position:
    Ny_E = [0 1 0]'; % selected y-axis direction
end

% Find x-axis of N (North):
Nx_E = cross(Ny_E,Nz_E); % Final axis found by right hand rule

% Form R_EN from the unit vectors:
R_EN = R_Ee'*[Nx_E Ny_E Nz_E]; % R_Ee selects correct E-axes, see R_Ee.m for details
