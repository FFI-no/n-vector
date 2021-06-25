function R_AB = xyz2R(x,y,z)
% xyz2R Creates a rotation matrix from 3 angles about new axes in the xyz order. 
%   R_AB = xyz2R(x,y,z) 
%   The rotation matrix R_AB is created based on 3 angles x,y,z about new
%   axes (intrinsic) in the order x-y-z. The angles (called Euler angles or
%   Tait-Bryan angles) are defined by the following procedure of successive
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
%   x,y,z [rad]	    Angles of rotation about new axes.
%
%   OUT:
%   R_AB  [no unit]	3x3 rotation matrix (direction cosine matrix) such that the 
%                   relation between a vector v decomposed in A and B is 
%                   given by: v_A = R_AB * v_B
% 
%   See also R2xyz, zyx2R, R2zyx.


%   This file is part of NavLab and is available from www.navlab.net/nvector
%
%   Copyright (c) 2021, Norwegian Defence Research Establishment (FFI)
%   All rights reserved.

% Originated: 1996.10.01 Kenneth Gade, FFI
% Modified:

cz = cos(z); sz = sin(z);
cy = cos(y); sy = sin(y);
cx = cos(x); sx = sin(x);

R_AB = [cy*cz              -cy*sz      sy
      sy*sx*cz+cx*sz -sy*sx*sz+cx*cz -cy*sx
      -sy*cx*cz+sx*sz sy*cx*sz+sx*cz  cy*cx];
  
