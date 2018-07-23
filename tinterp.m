% interpolate u1, u2 directly
% standard interpolation

function ts = tinterp(t, u1, u2)

ts = u1*(1 - t) + u2*t;

end