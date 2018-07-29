%% eigenvalues of operator

load chen12_1500;

% % double pulse
index = 2;
% y2 = ydF(:,index);
y2 = ydFD(:,index);

% default parameters
N = length(x);                  % current grid size
L = ceil(abs(x(1)));            % current domain length;
h = abs(x(2) - x(1));           % current grid spacing
c = par.c;
y = ynew;

% compute A0
[fval, A0] = ChenExp(x,y,par,D,config);

% Melnikov integral
Dy = D(:,:,1)*y;
M = sum( Dy.*Dy )*h;

% % Fourier
% D = DF;
% config = configF;
% x = xF;
% y = yF;
% h = x(2) - x(1);

Dy = D(:,:,1)*y;
M = sum( Dy.*Dy )*h;

%% eigenvalues

% compute approximate eigenvalues using formula from my research

% find peak distance
% to do this, we need a finer grid
Nfine = 3000;
Lfine = 30;

% % Finite difference
% xfine = linspace(-Lfine, Lfine, Nfine)';
% hfine = xfine(2) - xfine(1);
% y2fine = spline(x,y2,xfine);
% y2fine(isnan(y2fine)) = 0;
% yfine = spline(x,y,xfine);
% yfine(isnan(yfine)) = 0;

% [xfine, y2fine] = fsolve( @(y) ChenExp(xfine,yfine,par,Dfine,config), y, options);

% [pks, locs] = findpeaks(-y2fine);
% indices = find(abs(pks) > 1);
% peaklocs = locs(indices);
% peakdist = peaklocs(2) - peaklocs(1);
% L1 = peakdist * hfine / 2;

[pks, locs] = findpeaks(-y2);
indices = find(abs(pks) > 1);
peaklocs = locs(indices);
peakdist = peaklocs(2) - peaklocs(1);
L1 = peakdist * h / 2;

% get x position associated with this
% xpos = find(abs(xfine - L1) < 0.01);
% xneg = find(abs(xfine + L1) < 0.01);

xpos = find(abs(x - L1) < 0.01);
xneg = find(abs(x + L1) < 0.01);

% % q and its derivatives 
% D1q = spline(x,D(:,:,1)*y,xfine);
% D2q = spline(x,D(:,:,2)*y,xfine);
% D3q = spline(x,D(:,:,3)*y,xfine);
% D4q = spline(x,D(:,:,4)*y,xfine);

D1q = D(:,:,1)*y;
D2q = D(:,:,2)*y;
D3q = D(:,:,3)*y;
D4q = D(:,:,4)*y;

Qprime = [ D1q D2q D3q D4q ];
Psi = [ (-D4q - c^2 * D2q) (D3q + c^2 * D1q) -D2q D1q];

a = Qprime(xpos,:)*Psi(xneg,:)';

lest = sqrt(2*a / M);
lest0 = 2*a / M;

[fval, A0] = ChenExp(x,y2,par,D,config);
lambda0 = eig(A0);

