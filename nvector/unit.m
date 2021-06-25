function unit_vector = unit(vector)
%unit Makes input vector unit length, i.e. norm == 1.
%   unit_vector = unit(vector) 
%   makes the general 1xm vector a unit vector (norm == 1). A matrix of
%   vectors (nxm) is accepted as input.


%   This file is part of NavLab and is available from www.navlab.net/nvector
%
%   Copyright (c) 2021, Norwegian Defence Research Establishment (FFI)
%   All rights reserved.

% Originated: 2002.07.04 Kenneth Gade, FFI
% Modified:   2012.04.25 Kristian Svartveit, FFI: Vectorization
% Modified:   2015.02.19 Kristian Svartveit, FFI: Bugfix and speedup
% Modified:   2015.11.04 Kristian Svartveit, FFI: Fast for both matrix of
%                        vectors and single vector input

[m, n] = size(vector);

if n == 1 % Quicker solution for single vector input:
    % Find vector norm
    current_norm = norm(vector);
    if current_norm == 0 % If the vector has norm == 0, i.e. all elements in the
        % vector are zero, the unit vector [1 0 0]' is returned.
        unit_vector = zeros(m,1);
        unit_vector(1) = 1;
    else
        % Divide vector with the norm
        unit_vector = vector/current_norm;
    end
else % Quicker solution for vectorized input:
    % Find each vector norm
    current_norm = sqrt(sum(abs(vector).^2));
    
    % Divide each vector with its corresponding norm
    unit_vector = bsxfun(@rdivide, vector, current_norm);
    
    % If a vector has norm == 0, i.e. all elements in the vector are zero, the
    % unit vector [1 0 0]' is returned.
    zero_idx = current_norm == 0;
    single_unit_vector = zeros(m,1);
    single_unit_vector(1) = 1;
    unit_vector(:,zero_idx) = repmat(single_unit_vector,1,sum(zero_idx));
end
