% exact soliton solution to KdV
% from Pego and Weinstein

function u = KdVsoliton(x, par)

c = par.c;
p = par.p;

m = ( (1/2) * c * (p+1)*(p+2) )^(1/p);
x1 = (1/2)*p*sqrt(c)*x;

u = m * ( sech(x1) ).^(2/p);

end