function [x,y,z] = R2xyz(R_AB)
% R2xyz 3 angles about new axes in the xyz order are found from a rotation matrix. 
%   [x,y,z] = R2xyz(R_AB) 
%   3 angles x,y,z about new axes (intrinsic) in the order x-y-z are found
%   from the rotation matrix R_AB. The angles (called Euler angles or
%   Tait–Bryan angles) are defined by the following procedure of successive
%   rotations:
%   Given two arbitrary coordinate frames A and B. Consider a temporary frame
%   T that initially coincides with A. In order to make T align with B, we
%   first rotate T an angle x about its x-axis (common axis for both A and T).
%   Secondly, T is rotated an angle y about the NEW y-axis of T. Finally, T
%   is rotated an angle z about its NEWEST z-axis. The final orientation of
%   T now coincides with the orientation of B.
%
%   The signs of the angles are given by the directions of the axes and the
%   right hand rule.
%
%   IN:
%   R_AB  [no unit]	3x3 rotation matrix (direction cosine matrix) such that the
%                   relation between a vector v decomposed in A and B is
%                   given by: v_A = R_AB * v_B
%
%   OUT:
%   x,y,z [rad]	    Angles of rotation about new axes.
%
%   See also xyz2R, R2zyx, zyx2R.


%   This file is part of NavLab and is available from www.navlab.net/nvector
%
%   Copyright (c) 2021, Norwegian Defence Research Establishment (FFI)
%   All rights reserved.

% Originated: 1996.10.01 Kenneth Gade, FFI
% Modified: 2021.05.04 Kenneth Gade, FFI. Included handling of singularities


% cos_y is based on as many elements as possible, to average out
% numerical errors. It is selected as the positive square root since
% y: [-pi/2 pi/2]
cos_y = sqrt((R_AB(1,1)^2 + R_AB(1,2)^2 + R_AB(2,3)^2 + R_AB(3,3)^2)/2);

n_of_eps_to_define_singularity = 10;
% Check if (close to) zyx Euler angle singularity:
if cos_y > n_of_eps_to_define_singularity*eps
    
    % Outside singularity:        
    % atan2: [-pi pi]
    z = atan2(-R_AB(1,2),R_AB(1,1));
    x = atan2(-R_AB(2,3),R_AB(3,3));
    
    sin_y = R_AB(1,3);
    
    y = atan2(sin_y,cos_y);
    
else
    
    % In singularity (or close to), i.e. y = +pi/2 or -pi/2:
    y = sign(R_AB(1,3))*pi/2; % Selecting y = +-pi/2, with correct sign
    
    % Only the sum/difference of x and z is now given, choosing x = 0:
    x = 0;
    
    % Lower left 2x2 elements of R_AB now only consists of sin_z and cos_z.
    % Using the two whose signs are the same for both singularities:
    z = atan2(R_AB(2,1),R_AB(2,2));
    
end
