function [xout, uout, fval] = fsolveequation(eq, xold, uold, par, N, L, config, opts)

%% setup

% default options
options=optimset('Display','iter','Jacobian','on','MaxIter',100);
options.TolFun = 1e-12;
options.TolX = 1e-12;

% set options
if exist('opts','var')
    if isfield(opts, 'iter')
        options.MaxIter = opts.iter;
    end
    if isfield(opts, 'Jacobian')
        options.Jacobian = opts.Jacobian;
    end
end

% extract old versions of N and L 
N_old = length(xold);
L_old = ceil(abs(xold(1)));

% differentiation matrices
if strcmp(config.method,'Fourier')
    D = D_fourier(N, L, config.degree);
    
elseif strcmp(config.method,'Chebyshev')
    N_old = N_old + config.num_Dirichlet;
    [D, xout] = D_cheb(N, L, config.degree, config);
end

% if N is different from size of grid xold, 
% or L is different from domain of xold,
% then interpolate onto a new grid
if (N ~= N_old) || (L ~= L_old)
    % if we have a periodic domain, recompute xout
    % for Chebyshev, this was already done
    if strcmp(config.BC,'periodic')
        % for periodic domain, have to remove final point
        xout = linspace(-L, L, N+1)';
        xout = xout(1:end-1);
    end
    
    % interpolate onto new grid
    ustart = spline(xold,uold(1:end-1),xout);
    % if we have some NaN values at the ends of u, replace these with 0
    ustart(isnan(ustart)) = 0;

% otherwise keep the same grid for x and the same u
else
    xout = xold;
    ustart = uold;
end

%% solve nonlinear problem using fsolve

% call fsolve with these options

[uout,fval] = fsolve(@(u) eq(xout, u, par, D, config), ustart, options);


