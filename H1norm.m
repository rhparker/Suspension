% discrete H1 norm

function n = H1norm(x, y, D)

D1 = D(:,:,1);
h = x(2) - x(1);
integrand = y.^2 + (D1*y).^2;
s = h*sum(integrand);

n = s.^(1/2);

end