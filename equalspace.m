% make equally spaced in H1 norm
% naive method

function Unew = equalspace(x,U,N,D)
    frac = 0.1;
    initsp = H1spread(x,U,D);
    Unew = U;
    
    for outerindex = [1:10]
        for index = [2:N-1]
            sp = H1spread(x,Unew,D);
            Uleft = Unew;
            Uleft(:,index) = (1 - frac)*Unew(:,index) + frac*Unew(:,index-1);
            spleft = H1spread(x,Uleft,D);
            Uright = Unew;
            Uright(:,index) = (1 - frac)*Unew(:,index) + frac*Unew(:,index+1);
            spright = H1spread(x,Uright,D);
            m = min([sp spleft spright]);
            if (m == spleft)
%                 disp('left');
                Unew = Uleft;
            elseif (m == spright)
%                 disp('right');
                Unew = Uright;
            end
        end
    end

end