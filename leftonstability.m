
% grid for x
L = 20;
N = 512;

% starting parameters
par.c = 0;

% method configuration
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
% % config.Neumann = 'LR';
% config.num_Dirichlet = 2;
% [D, xout] = D_cheb(N, L, config.degree, config);

par.b = -2;

% generate c = 0 solution
u = lefton(xout, 0, 20, par);
fval = Bfamily(xout, u, par, D, config);

% disp(['norm when plugged in: ', num2str(norm(fval))]); 

% load solution
load uc_out;
index = 1;
% u = uc(1:end-1, index);
% par.c = uc(end, index);

% figure;
% plot(xout, fval);


%%

% opts.Jacobian = 'on';
% [~, uout, fval] = fsolveequation(@Bfamily_int, xout, u, par, N, L, config, opts);
% plot(xout, uout);

% opts.Jacobian = 'on';
% [~, uout, fval] = fsolveequation(@Bfamily, xout, uout, par, N, L, config, opts);

uout = u;

Id = eye(length(xout));
B = inv(Id - D(:,:,2));

[fval, J, Jnl] = Bfamily(xout,u,par,D,config);
    
M1 = par.c * D(:,:,1) + B*Jnl;

[V, lambdaD] = eig(M1);

lambda = diag(lambdaD);

figure;
plot(lambda, '.');
title(['Lefton Eigenvalues (Fourier spectral methods), c = 0']);
xlabel('real part');
ylabel('imaginary part');
