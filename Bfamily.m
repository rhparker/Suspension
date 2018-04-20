% B family of equations (equilibrium, comoving frame)

% x: x values (not needed here)
% u: value of (potential) solution
% c: parameters (b, c)
% D: differentiation matrices
% config: configuration

function [F,J,Jnl] = Bfamily(x,u,par,D,config)
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

% check to see if we passed the derivatives separately
if isfield(par, 'D1u')
   D1u = par.D1u; 
else
   D1u = D1*u;
end

if isfield(par, 'D2u')
   D2u = par.D2u; 
else
   D2u = D2*u;
end

if isfield(par, 'D3u')
   D3u = par.D3u; 
else
   D3u = D3*u;
end

% equation of B family

% linear operator for linear part
LN = -c*D1 + c*D3;

LNu = -c*D1u + c*D3u;

% linear and nonlinear part
F  = LNu + (b+1).*(u.*(D1u)) - u.*(D3u) - b.*(D1u).*(D2u);
%% Jacobian
if nargout > 1     

% % nonlinear terms
% NL1 = sparse(1:N,[1:N],D1*u,N,N) + sparse(1:N,[1:N],u,N,N)*D1;
% NL2 = sparse(1:N,[1:N],D3*u,N,N) + sparse(1:N,[1:N],u,N,N)*D3;
% NL3 = sparse(1:N,[1:N],D2*u,N,N)*D1 + sparse(1:N,[1:N],D1*u,N,N)*D2;
% 
% J = LN + (b+1).*NL1 - NL2 - b.*NL3;

U0   = sparse(1:N,[1:N],u,N,N);
DU0  = sparse(1:N,[1:N],D1u,N,N);
D2U0 = sparse(1:N,[1:N],D2u,N,N);
D3U0 = sparse(1:N,[1:N],D3u,N,N);

NL0 = -(b+1).*DU0 + D3U0;
NL1 = ( -(b+1).*U0 + b*D2U0 )*D1;
NL2 = b*DU0*D2;
NL3 = U0*D3;

Jnl = NL0 + NL1 + NL2 + NL3;
J = -LN + Jnl;

end
