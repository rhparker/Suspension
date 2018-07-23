% interpolate u1, u2 via s 
% as in (2.4) in Chamard (2011)

function us = sinterp(t, s, u1, u2)

if t <= 1/2
    us = s*(2*t) + u1*(1 - 2*t);
else
    us = u2*(2*t-1) + s*(2 - 2*t);
end

end