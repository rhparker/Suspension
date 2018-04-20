%% Centered finite difference differentiation matrices
% N is number of grid points
% h is mesh size
% num is number of derivatives
% BC is boundary conditions ('Neumann', 'periodic', 'Dirichlet')

function D = D_fdiff(N, h, num, BC, pts)

% default is 3-point stencil so we can get to 5th derivative
if ~exist('pts','var')
    pts = 3;
end

D = zeros([N, N, num]);

% compute differentiation matrices
for index = 1:num
    
    % specify stencil size
    if pts > 0
        M = fddiff(N, index, h, pts, BC);
    % take minimum stencil
    else
        M = fddiffeasy(N, index, h, BC);
    end

    D(:,:,index) = M;
end

end