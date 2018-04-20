%% generates finite difference matrix using Fornberg weights
%   N       : size of matrix
%   m       : order of derivative
%   h       : spatial step size
%   pts     : num of points to use on either side of center
%   BC      : boundary conditions ('periodic' for periodic)
% note that pts must be at least m/2 (rounded up) for this to work

function D = fddiff(N, m, h, pts, BC)
    w   = Fornberg_weights(0,-pts:pts,m);
    wts = w(m+1,:);
    D   = sparse_multidiag(N, wts, BC);
    stencil = 2*pts + 1;
    
    if (strcmp(BC,'Neumann'))
        % deal with ghost points
        D = D + ghost_points(N, wts);
        % if odd derivative, then first and last rows are 0
        % since odd derivatives are 0 at boundary
        if mod(m, 2) == 1
            D(1,:) = 0;
            D(N,:) = 0;
        end
        
     % Dirichlet BCs the easy way, zeroing out columns
     elseif (strcmp(BC,'Dirichlet'))
        % if even derivative, then first and last rows are 0
        % since even derivatives are 0 at boundary
        if mod(m, 2) == 0
            D(1,:) = 0;
            D(N,:) = 0;
        end
        % first and last grid points are always 0 (since Dirichlet)
        % so we zero out those columns
        D(:,1) = 0;
        D(:,N) = 0;
        
    % one-sided differences on the ends, no BCs enforced
    elseif (strcmp(BC,'OneSided'))
        noncentered = (stencil - 1)/2;
        % put noncentered differences on appropriate rows
        allweights = ufdwt(h,stencil,m);
        for index = 1:noncentered
            D(index, 1:stencil) = allweights(index, :);
            D(N - index + 1, end - stencil + 1:end) = allweights(stencil - index + 1, :);
        end
        
    % Dirichlet BCs using one-sided differences on the ends
    elseif (strcmp(BC,'DirichletOneSided'))
        noncentered = (stencil - 1)/2;
        % put noncentered differences on appropriate rows
        allweights = ufdwt(h,stencil,m);
        for index = 1:noncentered
            D(index, 1:stencil) = allweights(index, :);
            D(N - index + 1, end - stencil + 1:end) = allweights(stencil - index + 1, :);
        end
        % if even derivative, then first and last rows are 0
        % since even derivatives are 0 at boundary
        if mod(m, 2) == 0
            D(1,:) = 0;
            
            D(N,:) = 0;
        end
        % first and last grid points are always 0 (ince Dirichlet)
        % so we zero out those columns
        D(:,1) = 0;
        D(:,N) = 0;
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

% generates sparse matrix to add to diff matrix 
% for ghost points
%   N     : size of matrix
%   wts   : weights for multi-diagonal of matrix
% wts needs to have odd length
function S = ghost_points(N, wts)
    % start with zero sparse matrix
    S = sparse(N, N);
    % number of rows we need to deal with is number of 
    % weights which protrude outside of matrix
    % i.e. number of points on either side of center
    pts = (length(wts) - 1) / 2;
    % values at ghost points are added to those of corresponding
    % points flipped over 0 or N
    for i = 1:pts
        index = pts - i + 1;
        S(i, 2:index+1) = flip( wts(1:index) );
        S(N-i+1,N-index:N-1) = flip( wts( end-index+1:end) );
    end
end