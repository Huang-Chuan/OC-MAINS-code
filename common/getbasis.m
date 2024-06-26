function basis = getbasis(order)
    switch order
        case 1
            basis = [0 0 1;...
                     0 1 0;...
                     1 0 0;...
                    ];            
        case 2
            basis = [0 0 1;...
                     0 1 0;...
                     1 0 0;...
                     1 1 1;...
                     0 0 0];
        case 3
            basis =  [0 0 1;...
                     0 1 0;...
                     1 0 0;...
                     1 1 1;...
                     0 0 0;...
                     1 1 0;...
                     1 0 1;...
                    -1 0 1];
        otherwise
            disp('Not supported!')
    end
end