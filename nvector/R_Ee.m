function R_Ee_selected = R_Ee()
% R_Ee Selects axes of the coordinate frame E.
%   This file controls the axes of the coordinate frame E (Earth-Centred, 
%   Earth-Fixed, ECEF) used by the other files in this library
%
%   There are two choices of E-axes that are described in Table 2 in Gade
%   (2010):
%
%  * e: z-axis points to the North Pole and x-axis points to the point where
%       latitude = longitude = 0. This choice is very common in many fields.
%
%  * E: x-axis points to the North Pole, y-axis points towards longitude +90deg 
%       (east) and latitude = 0. This choice of axis directions ensures
%       that at zero latitude and longitude, N (North-East-Down) has the
%       same orientation as E. If roll/pitch/yaw are zero, also B (Body,
%       forward, starboard, down) has this orientation. In this manner, the
%       axes of E is chosen to correspond with the axes of N and B.
%
%   Based on this we get:
%   R_Ee = [0 0 1
%           0 1 0
%          -1 0 0]
%
%   The above R_Ee should be returned from this file when using z-axis to the 
%   North pole (which is most common). When using x-axis to the North 
%   pole, R_Ee should be set to I (identity matrix) (since the files in 
%   this library are originally written for this option).
%
%   Reference:
%   K Gade (2010): A Nonsingular Horizontal Position Representation, 
%   The Journal of Navigation, Volume 63, Issue 03, pp 395-417, July 2010. 
%   (www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)


%   This file is part of NavLab and is available from www.navlab.net/nvector
%
%   Copyright (c) 2021, Norwegian Defence Research Establishment (FFI)
%   All rights reserved.

% Originated: 2015.02.19 Kenneth Gade and Kristian Svartveit, FFI
% Modified:


% Select axes of E (by commenting/uncommenting the two choices):

% z-axis to the North Pole, default and most common choice:

R_Ee_selected = [0 0 1
                 0 1 0
                -1 0 0];

          
% x-axis to the North Pole, less common, but corresponds to typical choices of N and B:

% R_Ee_selected = [1 0 0
%                  0 1 0
%                  0 0 1];



