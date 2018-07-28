% essential spectrum gap

cmin = 0;
cmax = sqrt(2);
Nc = 200;
c = linspace(cmin, cmax, Nc);

rmin = -2;
rmax = 2;
Nr = 100;
r = linspace(rmin, rmax, Nr+1);

% plot of ess spec for c = 1
figure;
hold on;
lplus = r + sqrt(1 + r.^4);
lminus = r - sqrt(1 + r.^4);
plot(r, lplus, r, lminus);
plot(r, zeros(size(r)), '--');
legend('lambda+', 'lambda-', 'lambda = 0 line');
xlabel('r');
ylabel('imag part of (purely imaginary) lambda');

% find the minima for all our c
sqrterm = sqrt(1 + r.^4);
mins = zeros(size(c));
for index = [1 : length(c) ]
    mins(index) = min( c(index)*r + sqrterm );
end

figure;
plot(c,mins,'Linewidth',2);
xlabel('c');
ylabel('minimum imag part of lambda+');



