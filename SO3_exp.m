function R = SO3_exp(v)
    theta = norm(v);
    if theta < eps
        R = eye(3);
    else
        k = v / theta;
        K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
        R = eye(3) + sin(theta) * K + (1 - cos(theta)) * K^2;
    end
end
