
% grid for x
L = 20;
N = 256;

% starting parameters
par.c = 2;
par.b = -5;

% % method configuration
%  Fourier
config.method = 'Fourier';
config.BC = 'periodic';
config.degree = 3;
xout = linspace(-L, L, N+1)';
xout = xout(1:end-1);
D = D_fourier(N, L, config.degree);

% % Chebyshev
% config.method = 'Chebyshev';
% config.degree = 3;
% config.Dirichlet = 'LR';
% config.num_Dirichlet = 2;
% [D, xout] = D_cheb(N, L, config.degree, config);

u = lefton(xout, 0, 1, par);

opts.Jacobian = 'on';
[~, uout, fval] = fsolveequation(@Bfamily_int, xout, u, par, N, L, config, opts);
plot(xout, uout);

% opts.Jacobian = 'off';
% [~, uout, fval] = fsolveequation(@Bfamily, xout, uout, par, N, L, config, opts);

[fval, J] = Bfamily(xout,uout,par,D,config);
Id = eye(N);
B = inv(Id - D(:,:,2));
[V, lambdaD] = eig(B*J);
lambda = diag(lambdaD);
realL = lambda( find(abs(imag(lambda)) == 0) );
posL = realL( find(realL >= 0));