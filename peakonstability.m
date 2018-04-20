clear config;
clear par;

% grid for x
L = 15;
N = 512;

% starting parameters
par.c = 2;
par.b = 3;

% which method to use
% config.method = 'Fourier';
config.method = 'Chebyshev';

% method configuration

% for all methods
config.degree = 3;

% Fourier
if strcmp(config.method, 'Fourier')
    config.BC = 'periodic';
    config.num_Dirichlet = 0;
    xout = linspace(-L, L, N+1)';
    xout = xout(1:end-1);
    D = D_fourier(N, L, config.degree);

% Chebyshev
elseif strcmp(config.method, 'Chebyshev')
	% want odd number of grid points
    N = N + 1;
    config.Dirichlet = 'LR';
    config.num_Dirichlet = 2;
    [D, xout] = D_cheb(N, L, config.degree, config);
end
    
u = par.c.*exp(-abs(xout));

% options = optimset('Display','iter','Algorithm','levenberg-marquardt','MaxIter',500,'Jacobian','on');
% options.TolFun = 1e-12;
% options.TolX = 1e-12;
% [u1,fval,exitflag,output,jacobian1]  = fsolve( @(u) Bfamily_int(xout,u,par,D,config), u, options);


% % options to pass to fsolve
% opts.Jacobian = 'on';
% opts.iter = 1000;
% [~, uout, fval] = fsolveequation(@Bfamily_int, xout, u, par, N, L, config, opts);
% 
% plot(xout, uout,'.');
% plot(xout, u, xout, uout);
% 
% config.Neumann = 'R';
% [~, uout, fval] = fsolveequation(@Bfamily, xout, u, par, N, L, config, opts);
% plot(xout, u, xout, uout);

[fval, J] = Bfamily(xout,u,par,D,config);

Id = eye(N - config.num_Dirichlet);
B = inv(Id - D(:,:,2));
[V, lambdaD] = eig(B*J);
lambda = diag(lambdaD);
realL = lambda( find(abs(imag(lambda)) == 0) );
posL = realL( find(realL >= 0));

