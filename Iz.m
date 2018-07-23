% Functional L needed for the mountain pass algorithm

function result = Iz(x, z, par, D)

D1 = D(:,:,1);
D2 = D(:,:,2);
h = x(2) - x(1);
c = par.c;

% first, we need F, which is cumulative integral of 
% \tilde{f} = f(z + 1) = e^{z} - 1 
% from 0 to z
% we can just compute this!
% so we have F(z) = e^z - z - 1

% integrand
F = exp(z) - z - 1;
H = (1/2) * ( abs(D2*z).^2 - c^2 * abs(D1*z).^2 ) + F;
% H = abs(D2*z).^2 - c^2 * abs(D1*z).^2;

% integrate over x
result = sum( H(1:end-1) ) * h;

end