function R_AB = zyx2R(z,y,x)
% zyx2R Creates a rotation matrix from 3 angles about new axes in the zyx order. 
%   R_AB = zyx2R(z,y,x) 
%   The rotation matrix R_AB is created based on 3 angles z,y,x about new
%   axes (intrinsic) in the order z-y-x. The angles (called Euler angles or
%   Tait-Bryan angles) are defined by the following procedure of successive
%   rotations:
%   Given two arbitrary coordinate frames A and B. Consider a temporary frame 
%   T that initially coincides with A. In order to make T align with B, we 
%   first rotate T an angle z about its z-axis (common axis for both A and T). 
%   Secondly, T is rotated an angle y about the NEW y-axis of T. Finally, T 
%   is rotated an angle x about its NEWEST x-axis. The final orientation of 
%   T now coincides with the orientation of B.
%
%   The signs of the angles are given by the directions of the axes and the 
%   right hand rule.
%
%   Note that if A is a north-east-down frame and B is a body frame, we 
%   have that z = yaw, y = pitch and x = roll. 
%
%   IN: 
%   z,y,x [rad]	    Angles of rotation about new axes.
%
%   OUT:
%   R_AB  [no unit]	3x3 rotation matrix (direction cosine matrix) such that the 
%                   relation between a vector v decomposed in A and B is 
%                   given by: v_A = R_AB * v_B
% 
%   See also R2zyx, xyz2R, R2xyz.


%   This file is part of NavLab and is available from www.navlab.net/nvector
%
%   Copyright (c) 2021, Norwegian Defence Research Establishment (FFI)
%   All rights reserved.

% Originated: 1996.10.01 Kenneth Gade, FFI
% Modified:

cz = cos(z); sz = sin(z);
cy = cos(y); sy = sin(y);
cx = cos(x); sx = sin(x);

R_AB = [cz*cy -sz*cx+cz*sy*sx  sz*sx+cz*sy*cx  
        sz*cy  cz*cx+sz*sy*sx -cz*sx+sz*sy*cx  
         -sy       cy*sx            cy*cx    ];
        
