function dx=fp(t,x,p)
    switch p
        case 1
            dx = [x(1)/50*(-6*x(1)^2 - 6*x(2)^2 + 50);x(2)];
        case 2 
            dx = [-10*x(1) + 30/25*x(2)^3 - 15*x(2);-3*x(2)];    
    end
end