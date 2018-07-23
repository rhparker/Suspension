function p = betterpath(x, z1, z2, par, D, L)

% max of z1
Im = Iz(x, z1, par, D, L);
newIz = Im;

steps = 0;
max_steps = 50;
normdiff = norm( z2 - z1 );
threshold = 1e-5;
gOld = z2;

while (steps < max_steps) & (normdiff > threshold) & (newIz <= Im)
    gNew = (1/2)*(z1 + gOld);
    newIz = Iz(x, gNew, par, D, L);
    steps = steps + 1;
    normdiff = normdiff / 2;
    gOld = gNew;
end

fprintf('Steps to better max: %d \n',steps);


end