%% generate double pulses

% load data

load chen1322;

y = ynew;
c = par.c;
N = length(x);
L = -x(1);

%% double pulse construction

% paramaters from linearization about 0
% find the spatial eigenvalues
nu = roots([1 0 c^2 0 1]);
decay = abs(real(nu(1)));
freq  = abs(imag(nu(1)));


% find center
center = (N+1)/2;

% half wave
xhalf = x( center : end );
yhalf = y( center : end );

% find peaks of half wave
[pks, locs] = findpeaks(yhalf);
spacing = round( mean(locs(2:end) - locs(1:end-1)) / 4);
start = locs(1);

% which min/max we want
minmax = 2; 

join_pt = start + (minmax - 1) * spacing;

% right half-wave
yd_half    = y( join_pt : center + join_pt - 1 );
yd_right   = flip( yd_half(1:end-1) );
yd         = [ yd_half ; yd_right ];


%% fsolve

% run joined pulse through Newton solver
options = optimset('Display','iter','Jacobian','on');
[yd_out, fval] = fsolve( @(y) ChenExp(x,y,par,D,config), yd, options);

% plot double wave before and after Newton solver

lw = 1.5;

figure('DefaultAxesFontSize',16);
plot(x, yd, x, yd_out, 'Linewidth',lw);
% axis(ax);
legend('Initial guess','Output of fsolve');

y2 = yd_out;

%% Fourier

%  setup
N = 256;
configF = config;
configF.method = 'Fourier';
configF.BC = 'periodic';
configF.degree = 1;
xF = linspace(-L, L, N+1)';
xF = xF(1:end-1);
DF = D_fourier(N, L, config.degree);

% interpolate onto new grid
y2F = spline(x,y2,xF);
y2F(isnan(y2F)) = 0;

% run through Newton solver
options = optimset('Display','iter','Jacobian','on');
[y2F_out, fval] = fsolve( @(y) ChenExp(xF,y,par,DF,configF), y2F, options);
y2F = y2F_out;

% run single pulse through Newton solver
options = optimset('Display','iter','Jacobian','on');
yF = spline(x,y,xF);
yF(isnan(yF)) = 0;
[yF, fval] = fsolve( @(y) ChenExp(xF,y,par,DF,configF), yF, options);


%% eigenvalues

% a = 0.1;

% % Finite diff with weight
% Da = Dexp(D,a);
% [fval, A0] = ChenExp(x,y2,par,Da,config);
% D1 = Da(:,:,1);

% Finite diff
[fval, A0] = ChenExp(x,y2,par,D,config);
D1 = D(:,:,1);

% % Finite diff single pulse
% % [fval, A0] = ChenExp(x,y,par,D,config);
% % D1 = D(:,:,1);

% % Fourier
% [fval, A0] = ChenExp(xF,y2F,par,DF,configF);
% D1 = DF(:,:,1);

% Fourier with weight
% DFa = Dexp(DF,a);
% [fval, A0] = ChenExp(xF,y2F,par,DFa,configF);
% D1 = DFa(:,:,1);

% % Fourier single pulse
% [fval, A0] = ChenExp(xF,yF,par,DF,configF);
% D1 = DF(:,:,1);


N = length(D1);
A1 = -2*c * D1;
A2 = eye(N);
% EVP is A0 + lambda A1 + lambda^2 A2
% [V,lambda] = polyeig(A0,A1,A2);
[V, lambda] = quadeig(A2, A1, A0);

[V0, lambda0] = eig(A0);
lambda0 = diag(lambda0);

% plot showing essential spectrum boundary
rmin = -2;
rmax = 2;
xbound = 0.05;
Nr = 100;
r = linspace(rmin, rmax, Nr+1);
sbound = min(c*r + sqrt(1 + r.^4));
figure;
hold on;
plot(lambda, '.', 'MarkerSize', 10);
plot([0, 0], [-sbound, sbound], 'x', 'MarkerSize',15);
axis([-xbound xbound -sbound*1.5 sbound*1.5])
legend('Spectrum', 'Ess spec bound');
