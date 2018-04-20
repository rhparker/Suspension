function u = lefton(x, x0, k, par)

b = par.b;

% % this version is from Lafortune slides
% g = -(b + 1)/2;
% A = 2*( k^(b+1) ) / (1-b);
% argument = g*(x - x0);
% u = ( A*(1-b)/2 )*( cosh(argument).^(par.b / g) );

% this version is from the Degasparis et al (2003) paper
argument = ( (b+1)/2 ) * (x - x0);
exponent = 2/(b+1);
u = cosh(argument).^exponent;

end