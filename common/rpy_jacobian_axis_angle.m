function J = rpy_jacobian_axis_angle(axisAngle)
    u = axisAngle(1:3);
    na = axisAngle(end);
    na3 = na^3;
    t = abs(axisAngle(end));
    a = u * t;
    % First-order approximation of Jacobian wrt u, t.
    Jr = [t/(t^2*u(1)^2 + 1), 0, 0, u(1)/(t^2*u(1)^2 + 1);
          0, t/sqrt(1 - t^2*u(2)^2), 0, u(2)/sqrt(1 - t^2*u(2)^2);
          0, 0, t/(t^2*u(3)^2 + 1), u(3)/(t^2*u(3)^2 + 1)];

    % Jacobian of u, t wrt a.
    Ja = [(a(2)^2 + a(3)^2)/na3,        -(a(1)*a(2))/na3,        -(a(1)*a(3))/na3;
         -(a(1)*a(2))/na3,         (a(1)^2 + a(3)^2)/na3,        -(a(2)*a(3))/na3; 
         -(a(1)*a(3))/na3,              -(a(2)*a(3))/na3,   (a(1)^2 + a(2)^2)/na3;    
                a(1)/na,                 a(2)/na,                 a(3)/na];

    J = Jr * Ja;
end