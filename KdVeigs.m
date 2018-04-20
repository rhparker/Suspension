clear config;
clear par;

% grid for x
L = 20;

% starting parameters
par.p = 1;
par.c = 1;

% for all methods
config.degree = 3;

% %  Fourier
% N = 512;
% config.method = 'Fourier';
% config.BC = 'periodic';
% config.degree = 3;
% xout = linspace(-L, L, N+1)';
% xout = xout(1:end-1);
% D = D_fourier(N, L, config.degree);
% u = KdVsoliton(xout, par);
% [fval, J] = KdV3(xout, u, par, D, config);
% [V, lambdaD] = eig(J);
% lambda = diag(lambdaD);
% plot(lambda, '.');
% 
% % test for even and odd
% threshold = 0.1;
% V = V(2:end,:);
% evens = [];
% odds = [];
% for index = 1:N
%     if norm( V(:,index) - flip(V(:,index)) ) < threshold
%         disp( ['even eigenfunction: ' num2str(index) ])
%         evens = [evens index];
%     end
%     if norm( V(:,index) + flip(V(:,index)) ) < threshold
%         disp( ['odd eigenfunction: ', num2str(index) ])
%         odds = [odds index];
%     end
% end


% % finite difference, periodic
% config.method = 'FD';
% BC = 'periodic';
% % differentiation matrices
% N = 512;
% xout = linspace(-L, L, N+1)';
% h = xout(2) - xout(1);
% xout = xout(1:end-1);
% D = D_fdiff(N, h, config.degree, BC, 3);
% u = KdVsoliton(xout, par);

% finite difference, nonperiodic
config.method = 'FD';

% boundary conditions
% BC = 'Neumann';
BC = 'Dirichlet';

% differentiation matrices
N = 1000;
N = N+1;

% full line
xout = linspace(-L, L, N)';
h = xout(2) - xout(1);
D = D_fdiff(N, h, config.degree, BC, 0);

% soliton solution
u = KdVsoliton(xout, par);
Du = D(:,:,1)*u;

% options = optimset('Display','iter','Algorithm','levenberg-marquardt','MaxIter',500,'Jacobian','on');
% options.TolFun = 1e-12;
% options.TolX = 1e-12;
% [u1,fval,exitflag,output,jacobian1]  = fsolve( @(u) KdV3(xout,u,par,D,config), u, options);


% solution and Jacobian for full line
[fval, J] = KdV3(xout, u, par, D, config);


%%

% half line
xhalf = xout( ((N-1)/2)+1 : end );
uhalf = u( ((N-1)/2)+1 : end );
Duhalf = Du( ((N-1)/2)+1 : end );
Nhalf = (N+1)/2;
h = xhalf(2) - xhalf(1);
D = D_fdiff(Nhalf, h, config.degree, BC, 0);

% for simplicity
D1 = D(:,:,1);
D2 = D(:,:,2);
D3 = D(:,:,3);

% solution and Jacobian for half line
par.Du = Duhalf;
[fval, Jhalf] = KdV3(xhalf, uhalf, par, D, config);

%%

% deal with Dirichlet BCs
if (strcmp(BC, 'Dirichlet'))
    Jhalf = Jhalf(2:end-1, 2:end-1);
end

% full line
[V, lambdaD] = eig(J);

% half line
% [V, lambdaD] = eig(Jhalf);

lambda = diag(lambdaD);
plot(lambda, '.');
