%% Chebyshev differentiation matrices
% this version has no boundary conditions
% for interval [-L, L]

function [D, x] = D_chebNoBC(N, L, num, config)

% generate Chebyshev matrices for num derivatives
[x, DM] = chebdif(N, num);

% scale for length
if isfield(config,'half_line')
    % scale factor for [0, L]
    scale = 2/L;
    x = (L/2).*(x+1);
else
    % scale factor for [-L, L]
    scale = 1/L;
    x = L.*x;
end

D = zeros([N, N, num]);

% scale matrices appropriately
for index = 1:num
    M = DM(:,:,index);
    % scale to interval [-L, L]
    M = M * (scale^index);
    % put into array D
    D(:,:,index) = M;
end

end