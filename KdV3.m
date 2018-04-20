% 3rd order KdV

function [F,J] = KdV3(x,u,par,D,config)
% returns the right-hand side of our equation

%% operator

% paramaters (for convenience)
p = par.p;
c = par.c;

% identity operator
N = length(x);
Id = eye(N);

% differentiation operators, for convenience
D1 = D(:,:,1);
D2 = D(:,:,2);
D3 = D(:,:,3);

% linear part
LN = D3 - c*D1;

% linear and nonlinear part
F  = LN*u + (u.^p).*(D1*u);

%% Jacobian
if nargout > 1     
    
% check to see if we passed the derivative separately
if isfield(par, 'Du')
   Du = par.Du; 
else
   Du = D1*u;
end

% nonlinear terms
NL1 = p * sparse(1:N,[1:N],u.^(p-1),N,N) * sparse(1:N,[1:N],Du,N,N);
NL2 = sparse(1:N,[1:N],u.^p,N,N) * D1;
J = LN + NL1 + NL2;
 
end
