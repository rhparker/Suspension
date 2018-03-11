% B family of equations (equilibrium, comoving frame)

% x: x values (not needed here)
% u: value of (potential) solution
% c: parameters (b, c)
% D: differentiation matrices
% config: configuration

function [F,J] = Bfamily(x,u,par,D,config)
% returns the right-hand side of our equation

%% operator

% paramaters (for convenience)
b = par.b;
c = par.c;

% grid size parameter
N = length(x);

% differentiation operators, for convenience
D1 = D(:,:,1);
D2 = D(:,:,2);
D3 = D(:,:,3);

% equation of B family

% linear operator for linear part
LN = -c*D1 + c*D3;

% linear and nonlinear part
F  = LN*u + (b+1).*(u.*(D1*u)) - u.*(D3*u) - b.*(D1*u).*(D2*u);
%% Jacobian
if nargout > 1     

% nonlinear terms
NL1 = sparse(1:N,[1:N],D1*u,N,N) + sparse(1:N,[1:N],u,N,N)*D1;
NL2 = sparse(1:N,[1:N],D3*u,N,N) + sparse(1:N,[1:N],u,N,N)*D3;
NL3 = sparse(1:N,[1:N],D2*u,N,N)*D1 + sparse(1:N,[1:N],D1*u,N,N)*D2;

J = LN + (b+1).*NL1 - NL2 - b.*NL3;
 
end
