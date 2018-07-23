% minimize energy

function zm = minenergy3(z, v, hmax, x, par, D)

steps = 20;
energy = Iz(x, z, par, D);
curr_energy = energy;
zm = z;

for index = flip([0:steps])
    znew = z + hmax * 2^(-index) * v;
    newenergy = Iz(x, znew, par, D);
    disp(['step size: 2^-' num2str(index) '   energy diff: ' num2str(newenergy - energy)]);
    if (newenergy < curr_energy)
        zm = znew;
        curr_energy = newenergy;
    else
        break;
    end
end

end