function deg_angle = deg(rad_angle)
%deg Converts angle in radians to degrees.
%   deg_angle = deg(rad_angle)
% 
%   IN:
%   angle in radians
%
%   OUT:
%   angle in degrees
%
%   See also rad.


%   This file is part of NavLab and is available from www.navlab.net/nvector
%
%   Copyright (c) 2021, Norwegian Defence Research Establishment (FFI)
%   All rights reserved.

% Originated: 1996 Kenneth Gade, FFI

deg_angle = rad_angle*180/pi;


