function r=PosMagArray()
% Output
%        r:  position of the sensors    --- 3 x (number of sensors) vector
    r = NaN(3, 30);
%   Comment out the below and replace r with the correponding positions 
%   if you want to use your own configuration
    dx = 0.05;
    dy = 0.05;
    kk = 0;

    for jj=1:5
        for ii=1:6
            kk=kk+1;
            if ii<4
                r(1, kk) = (ii-3.5)*dx-0.5*dx;
            else
                r(1, kk) = (ii-3.5)*dx+0.5*dx;
            end
            r(2, kk) = -(jj-3)*dy;
            r(3, kk) = 0;
        end
    end
end
