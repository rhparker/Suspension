clear config;

% grid for x
L = 1;
N = 8;

% starting parameters
par.c = 2;
par.b = 1;

% which method to use
config.method = 'Chebyshev';

% method configuration

% for all methods
config.degree = 3;

% want odd number of grid points
N = N + 1;
% config.Dirichlet = 'LR';
config.Neumann = 'LR';

% differentiation matrices
[D, xout] = D_chebNoBC(N, L, config.degree, config);
Id = eye(length(xout));
Z = Id;

% boundary conditions
BC = [];

% deal with Dirichlet BCs, if any
if isfield(config, 'Dirichlet')
    % Neumann on R
    if strcmp(config.Dirichlet,'R')
        BC = [BC;  Id(1,:) ];        
    elseif strcmp(config.Dirichlet,'L')
        BC = [BC ; Id(end, :) ];        
    elseif strcmp(config.Dirichlet,'LR')
        BC = [BC; Id(1,:); Id(end,:) ];        
    end
end

% deal with Neumann BCs, if any
if isfield(config, 'Neumann')
    % Neumann on R
    if strcmp(config.Neumann,'R')
        BC = [ BC; D(1,:,1) ];        
    elseif strcmp(config.Neumann,'L')
        BC = [ BC; D(end,:,1) ];        
    elseif strcmp(config.Neumann,'LR')
        BC = [ BC; D(1,:,1); D(end,:,1) ];        
    end
end

Z = null(BC);
B = Z' * D(:,:,1) * Z;

u = sin(pi * xout/ 2);

l = eig(B);
plot(l, '.');
