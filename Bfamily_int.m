% B family of equations (equilibrium, comoving frame)

% x: x values (not needed here)
% u: value of (potential) solution
% c: parameters (b, c)
% D: differentiation matrices
% config: configuration
% equation for B family, integrated onces

function [F,J] = Bfamily_int(x,u,par,D,config)
% returns the right-hand side of our equation

%% operator

% paramaters (for convenience)
b = par.b;
c = par.c;

% identity operator
N = length(x);
Id = eye(N);

% differentiation operators, for convenience
D1 = D(:,:,1);
D2 = D(:,:,2);

% equation of B family

% linear operator for linear part
LN = -c*Id + c*D2;

% linear and nonlinear part
F  = LN*u + ((b+1)/2).*(u.^2) + ((1-b)/2).*((D1*u).^2) - u.*(D2*u);

%% Jacobian
if nargout > 1     

% nonlinear terms
NL1 = (b+1).*sparse(1:N,[1:N],u,N,N) + (1-b).*sparse(1:N,[1:N],D1*u,N,N)*D1;
NL2 = -sparse(1:N,[1:N],D2*u,N,N) - sparse(1:N,[1:N],u,N,N)*D2;
J = LN + NL1 + NL2;
 
end
