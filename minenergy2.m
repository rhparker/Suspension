% minimize energy

function zm = minenergy2(z, v, hmax, x, par, D)

steps = 20;
energy = Iz(x, z, par, D);
disp(['starting energy: ' num2str(energy)]);
zm = z;

for index = [0:steps]
    znew = z + hmax * 0.5^index * v;
    newenergy = Iz(x, znew, par, D);
    disp(['step: ' num2str(index) 'energy diff: ' num2str(newenergy - energy)]);
    if (newenergy < energy)
        zm = znew;
        break;
    end
end

end