% finite difference, nonperiodic

config.method = 'FD';
% BC = 'Neumann';
BC = 'Dirichlet';

% differentiation matrices
L = 1;
N = 101;
x = linspace(0, L, N)';
h = x(2) - x(1);

s = sin(x);
c = cos(x);

D = D_fdiff(N, h, 3, BC, 0);

D1 = D(:,:,1);
D2 = D(:,:,2);
D3 = D(:,:,3);

% for Dirichlet
u = (x.^3) .* (x - L).^3;

% for Neumann
% u = (x.^4).*(x - L).^4 + 5;

test = D3*u;

plot(x, test);