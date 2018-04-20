clear config;

% grid for x
L = 20;         % domain length
N = 512;        % half line points

% want odd number of grid points
N = N + 1;

% full line
xfull = linspace(-L, L, 2*N - 1)';

% half line
xhalf = linspace(0, L, N)';

% grid spacing
h = xhalf(2) - xhalf(1);

% starting parameters
par.c = 2;
par.b = 2;

% which method to use
config.method = 'FD';

% for all methods, need three derivatives
config.degree = 3;

% use "exact solution"; do not fsolve
uhalf = par.c.*exp(-abs(xhalf));

% interpolate backwards from 0 so we can "differentiate"
offset = 10;
xinterp = xfull( N - offset + 1 : end );
uinterp = spline(xhalf,uhalf,xinterp);

% take "derivatives" of right half
D = D_fdiff(length(xinterp), h, config.degree, 'Neumann', 0);
D1u = D(:,:,1)*uinterp;
D2u = D(:,:,2)*uinterp;
D3u = D(:,:,3)*uinterp;
D1u = D1u(offset:end);
D2u = D2u(offset:end);
D3u = D3u(offset:end);

% convert these back to full derivatives
% odd derivatives odd, even derivatives even
D1full = [ -flip( D1u(2:end) ) ; D1u ];
D2full = [  flip( D2u(2:end) ) ; D2u ];
D3full = [ -flip( D3u(2:end) ) ; D3u ];

% pass these derivatives as parameters to the B-family
par.D1u = D1full;
par.D2u = D2full;
par.D3u = D3full;

% generate differentation operators for full line

% full line
x = xfull;

% exact solution  on x
u = (par.c).*exp(-abs(x));

% % BCs to use
BC = 'Dirichlet';
% BC = 'Neumann';

% differentiation matrices
D = D_fdiff( length(x), h, config.degree, BC, 0);

Id = eye(length(x));
B = inv(Id - D(:,:,2));

% use zero solution
% u = zeros(length(xout), 1);

%%

bvals = (2.5);
% bvals = linspace(0.5, 2.5, 9);

evals = zeros(length(bvals), length(x));

figure;
hold on;

for index = 1:length(bvals)
    % get Jacobian 
    par.b = bvals(index);
    
    [fval, J, Jnl] = Bfamily(x,u,par,D,config);
    
    M1 = par.c * D(:,:,1) + B*Jnl;

    [V, lambdaD] = eig(M1);
    lambda = diag(lambdaD);
    evals(index, :) = lambda;
    plot(lambda, '.');
    
%     realL = lambda( find(abs(imag(lambda)) == 0) );
%     posL = realL( find(realL >= 0));
end

%%

figure
hold on;
for index = 1:length(bvals)
    plot(evals(index, :), '.');
end
axis([-0.4, 0.5, -200, 200]);
title(['Peakon Eigenvalues, BCs: ' BC ]);
legendCell = cellstr(num2str(bvals', 'b=%.2f'))
legend(legendCell);
xlabel('real part');
ylabel('imaginary part');
