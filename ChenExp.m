% 3rd order KdV

function [F,J] = ChenExp(x,u,par,D,config)
% returns the right-hand side of our equation

%% operator

% paramaters (for convenience)
c = par.c;

% step size
h = x(2) - x(1);

% identity operator
N = length(x);
Id = eye(N);

% differentiation operators, for convenience
D1 = D(:,:,1);
D2 = D(:,:,2);
D3 = D(:,:,3);
D4 = D(:,:,4);

% linear part
LN = D4 + (c^2)*D2;

% linear and nonlinear part
% exponential version
F  = LN*u + ( exp(u) - 1 );

% Heaviside version
% F  = LN*u + ( max(u+1, 0) - 1 );

%% Jacobian
if nargout > 1     
    
% nonlinear terms
NL = sparse(1:N,[1:N],exp(u),N,N);
J = LN + NL;

% %% enforce norm, if requested
% % actually, enforce discrete L2 norm squared, since easier
% 
% if isfield(config, 'norm')
%     M = config.norm;
%     L2normSq = 1e-3*( sum( u.^2 ) - M ) ;
%     
%     % add element to function
%     F = [F; L2normSq];
%     
%     % add row/col to Jacobian
%     Jnormrow = 1e-3*2*u';
%     J = [J ; Jnormrow];
%     
% end
%  
end
