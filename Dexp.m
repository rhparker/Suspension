% takes diff operators and generates exp weight diff operators
function XD = Dexp(D, a)
    s = size(D);
    num = s(3);
    I = speye(s(1));
    
    XD = D;
    
    if (num >= 1)
        D1 = D(:,:,1);
        XD(:,:,1) = D1 - a*I;
    end
    if (num >= 2)
        D2 = D(:,:,2);
        XD(:,:,2) = D2 - 2*a*D1  + (a^2)*I;
    end
    if (num >= 3)
        D3 = D(:,:,3);
        XD(:,:,3) = D3 - 3*a*D2 + 3*(a^2)*D1  - (a^3)*I;
    end
    if (num >= 4)
        D4 = D(:,:,4);
        XD(:,:,4) = D4 - 4*a*D3 + 6*(a^2)*D2  - 4*(a^3)*D1 + (a^4)*I;
    end
    if (num >= 5)
        D5 = D(:,:,5);
        XD(:,:,5) = D5 - 5*a*D4 + 10*(a^2)*D3 - 10*(a^3)*D2 + 5*(a^4)*D1 - (a^5)*I;
    end
    
end