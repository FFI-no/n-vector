function rad_angle = rad(deg_angle)
%rad Converts angle in degrees to radians.
%   rad_angle = rad(deg_angle)
% 
%   IN:
%   angle in degrees
%
%   OUT:
%   angle in radians
%
%   See also deg.


%   This file is part of NavLab and is available from www.navlab.net/nvector
%
%   Copyright (c) 2021, Norwegian Defence Research Establishment (FFI)
%   All rights reserved.

% Originated: 1996 Kenneth Gade, FFI

rad_angle = deg_angle*pi/180;


