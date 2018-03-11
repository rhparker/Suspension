function u = lefton(x, x0, k, par)

b = par.b;
A = 2*( k^(par.b+1))/(1-par.b);
gamma = -(par.b + 1)/2;
argument = gamma*(x - x0);
u = ( A*(1-b)/2 )*( cosh(argument).^(par.b / gamma) );

end