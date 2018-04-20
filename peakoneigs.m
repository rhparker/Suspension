clear config;

% grid for x
L = 10;
N = 512;

% starting parameters
par.c = 2;
par.b = 3;

% which method to use
% config.method = 'Chebyshev';
config.method = 'FD';

% for all methods
config.degree = 3;

% want odd number of grid points
N = N + 1;

config.Dirichlet = 'LR';
config.Neumann = 'R';
config.half_line = 'True';

% differentiation matrices
% [D, xout] = D_chebNoBC(N, L, config.degree, config);
[D, xout] = D_cheb(N, L, config.degree, config);

Id = eye(length(xout));
B = Id - D(:,:,2);
Z = Id;

% boundary conditions
BC = [];
% deal with Dirichlet BCs, if any

% % deal with Dirichlet BCs, if any
% if isfield(config, 'Dirichlet')
%     % Neumann on R
%     if strcmp(config.Dirichlet,'R')
%         BC = [BC;  Id(1,:) ];        
%     elseif strcmp(config.Dirichlet,'L')
%         BC = [BC ; Id(end, :) ];        
%     elseif strcmp(config.Dirichlet,'LR')
%         BC = [BC; Id(1,:); Id(end,:) ];        
%     end
% end
% 
% % deal with Neumann BCs, if any
% if isfield(config, 'Neumann')
%     % Neumann on R
%     if strcmp(config.Neumann,'R')
%         BC = [BC; D(1,:,1) ];        
%     elseif strcmp(config.Neumann,'L')
%         BC = [BC; D(end,:,1) ];        
%     elseif strcmp(config.Neumann,'LR')
%         BC = [ BC; D(1,:,1); D(end,:,1) ];        
%     end
% end
% 
% Z = null(BC);
% B = inv( Z' * B * Z);

% use "exact solution"; do not fsolve
u = par.c.*exp(-abs(xout));

% plot(xout, u);

% we want to do this for many values of b
% bvals = linspace(1,3,11);

bvals = (2);

evals = zeros(length(bvals), length(B));

figure;
hold on;

for index = 1:length(bvals)
    % get Jacobian 
    par.b = bvals(index);
    [fval, J] = Bfamily(xout,u,par,D,config);

    [V, lambdaD] = eig(B*Z'*J*Z);
    lambda = diag(lambdaD);
    evals(index, :) = lambda;
    plot(lambda, '.');
    
%     realL = lambda( find(abs(imag(lambda)) == 0) );
%     posL = realL( find(realL >= 0));
end

