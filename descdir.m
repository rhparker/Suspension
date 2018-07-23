% steepest decent direction
% from Chen & McKenna (1997)

function v = descdir(z, K)

g = z - ftilde(z);
vbar = K\g;
v = (vbar - z)/norm(vbar - z);

end

% function f from the PDE
function f = ftilde(z)
    f = exp(z) - 1;
end