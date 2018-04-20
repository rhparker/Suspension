%% generates finite difference matrix using Fornberg weights
%   N       : size of matrix
%   m       : order of derivative
%   h       : spatial step size
%   BC      : boundary conditions ('periodic' for periodic)

function D = fddiffeasy(N, m, h, BC)
    % number of weights depends on order of derivative desired
    % we will use the simplest possible one
    pts = ceil(m/2);
    % Fornberg weights for centered FDs
    w   = Fornberg_weights(0,-pts:pts,m);
    wts = w(m+1,:);
    D   = sparse_multidiag(N, wts, BC);
    
    if (strcmp(BC,'Neumann'))
        % deal with ghost points
        if m == 2
            D(1,2) = D(1,2) + 1;
            D(end, end-1) = D(end, end-1) + 1;
        elseif m == 3
            D(2,2) = D(2,2) - 1/2;
            D(end-1, end-1) = D(end-1, end-1) + 1/2;
        end
 
        % if odd derivative, then first and last rows are 0
        % since odd derivatives are 0 at boundary
        if mod(m, 2) == 1
            D(1,:) = 0;
            D(N,:) = 0;
        end
        
     % Dirichlet BCs the easy way, zeroing out columns
     elseif (strcmp(BC,'Dirichlet'))
        % deal with ghost points
        if m == 1
            D(1,1:2) = D(1,1:2) + (-1/2)*[2 -1];
            D(end,end-1:end) = D(end,end-1:end) + (1/2)*[-1 2];
        elseif m == 3
            D(1,1:3) = D(1,1:3) + (-1/2)*[2 0 -1];
            D(end,end-2:end) = D(end,end-2:end) + (1/2)*[-1 0 2];
            D(1,1:2) = D(1,1:2) + 1*[2 -1];
            D(end,end-1:end) = D(end,end-1:end) + (-1)*[-1 2];
            D(2,1:2) = D(2,1:2) + (-1/2)*[2 -1];
            D(end-1,end-1:end) = D(end-1,end-1:end) + (1/2)*[-1 2];
        end
        % if even derivative, then first and last rows are 0
        % since even derivatives are 0 at boundary
        if mod(m, 2) == 0
            D(1,:) = 0;
            D(N,:) = 0;
        end
%         % first and last grid points are always 0 (since Dirichlet)
%         % so we zero out those columns
%         D(:,1) = 0;
%         D(:,N) = 0;
        
    end
    % divide by scaled spatial step size
    D = D / (h^m);
end

% generates sparse multidiagonal matrix S
%   N      : size of matrix
%   vals   : value to put on the diagonals
%   config : if 'periodic' then make periodic matrix
% vals must have an odd number of elements; the middle 
% element goes on the main diagonal
function S = sparse_multidiag(N, vals, config)
    len = length(vals);
    % start with zero sparse matrix
    S = sparse(N,N);
    % if odd number of vals, keep going
    if mod(len,2) == 1
        diags = -(len-1)/2 : (len-1)/2;
        for i = 1:len
            S = S + sparse_diag(N, diags(i), vals(i), config);
        end
    end
end

% generates sparse diagonal matrix S
%   N      : size of matrix
%   d      : which diagonal to use (0 is main diagonal)
%   val    : value to put on the diagonal
%   config : if 'periodic' then make periodic matrix
function S = sparse_diag(N, d, val, config)
    diag = abs(d);
    S = sparse(1:N-diag,diag+1:N,val*ones(N-diag,1),N,N);
    % if we want a periodic matrix
    if exist('config','var')
        if (strcmp(config,'periodic')) && (diag ~= 0)
            S = S + sparse_diag(N, N-diag, val)';
        end
    end
    % if we specify a negative diagonal, take transpose
    if d < 0
        S = S';
    end
end
