clear all;

% grid for x
L = 50;

% starting parameters
par.c = 1.354;
% par.c = 1.40;
par.p = 1;

% for convenience
c = par.c;

% for all methods
config.degree = 4;

% %  Fourier
% N = 513;
% config.method = 'Fourier';
% config.BC = 'periodic';
% config.degree = 4;
% x = linspace(-L, L, N+1)';
% x = x(1:end-1);
% D = D_fourier(N, L, config.degree);
% D2 = D(:,:,2);
% D4 = D(:,:,4);
% Id = eye(N);

% Finite difference
N = 513;
config.method = 'FD';
config.BC = 'none';
config.degree = 4;
x = linspace(-L, L, N)';
h = x(2) - x(1);
D = D_fdiff(N, h, config.degree, config.BC, 0);
Id = eye(N);
D1 = D(:,:,1);
D2 = D(:,:,2);
D4 = D(:,:,4);
D4(1,1) = 5 / h^4;
D4(end,end) = 5 / h^4;

% parameters
rho = acos(-c^2 / 2)/2;
tau = sin(rho);
sig = cos(rho);

% exact solution for Heaviside case
% from Chen & McKenna (1997)
% prob not necessary
z1 = exp(sig*x).*cos(tau*x);
z2 = (x.^2 / (2*c^2)) + cos(c*x);
% flip about x = 0
center = (N+1)/2;
z1 = [ z1(1:center) ; flip(z1(1:center-1)) ];
z2 = [ z2(1:center) ; flip(z2(1:center-1)) ];
z1 = -z1;
% plot(z1);

% this thing has I(ze) < 0
% anything should work, theoretically
% this follows Lemma 2.5 in McKenna & Chen
lambda = 0.5;
A = 8.5;

z = -exp(-(lambda * x).^2);
H = ( abs(D2*z).^2 - c^2 * abs(D1*z).^2 );
Mz1 = sum( H(1:end-1) ) * h;

z2 = -A * exp(-(lambda * x).^2);
Mze = Iz(x, z2, par, D);

% for the other one, just use the zero solution
z1 = zeros(size(z2));


%% String method (sort of)

% initialize
Nt = 10;
Ti = linspace(0, 1, Nt+1);
Ti = Ti(2:end);

U = tinterp(Ti,z1,z2);

grads = zeros(size(U));
energies = zeros(size(Ti));
hmax = 100;

% matrix for steepest descent direction
K = D4 + c^2*D2 + Id;

% % redo grid
% options = optimset('Display','iter','TolX', 1e0,'TolFun', 1e0);
% snew = fminsearch( @(s) H1spread(Ti,s,z1,z2,x,D), s, options);
% U = sinterp(x,snew,z1,z2);


%% iterate

maxreps = 10;

for rep = [1:maxreps]
    disp(['Iteration: ' num2str(rep)]);

    for index = [1:Nt]
        disp([' String: ' num2str(index)]);
        energies(index) = Iz(x, U(:,index), par, D);
        grads(:,index) = descdir(U(:,index),K);
        Unew(:,index) = minenergy3(U(:,index), grads(:,index), hmax, x, par, D);
    end
    
    U = Unew;
    disp('starting redist');
    U = equalspace(x,U,Nt,D);
    disp('redist complete');
    
end

%%

y0 = U(:,3);
options = optimset('Display','iter','Jacobian','on');
[ynew, fval] = fsolve( @(y) ChenExp(x,y,par,D,config), y0, options);


%% Mountain Pass algorithm
% 
% Nbar = 10;
% Ivals = zeros([Nbar,1]);
% 
% 
% for index = [1:Nbar]
%     Ivals(index) = Iz(x, ze * (index/Nbar), par, D, L);
% end
% 
% [M, iMax] = max(Ivals);
% 
% % try to find better path
% betterpath( x, ze * (iMax/Nbar), ze * ((iMax+1)/Nbar), par, D, L );



