% minimize energy

function zm = minenergy(z, v, hmax, x, par, D)

steps = 100;
t = linspace(0,hmax,steps+1);
energies = zeros(size(t));

% go in direction of gradient
vsteps = v*t;
zsteps = repmat(z,1,steps+1);
zsteps = zsteps + vsteps;

% compute energies
for index = [1:steps+1]
    energies(index) = Iz(x, zsteps(:,index), par, D);
end

[M, Ind] = min(energies);
plot(energies);
zm = zsteps(:,Ind);

end