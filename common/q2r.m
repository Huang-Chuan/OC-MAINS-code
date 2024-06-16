function rotm = q2r(q)
% Q2R Convert quaternion to rotation matrix
% input: q = [q_w, q_x, q_y, q_z]
% output: rotm = [3x3] rotation matrix

    % normalized the quaternion
    q = q/norm(q);

    q_w = q(1);
    q_x = q(2);
    q_y = q(3);
    q_z = q(4);


    rotm = [q_w^2+q_x^2-q_y^2-q_z^2,2*(q_x*q_y-q_w*q_z),2*(q_x*q_z+q_w*q_y);
           2*(q_x*q_y+q_w*q_z), q_w^2-q_x^2+q_y^2-q_z^2,2*(q_y*q_z-q_w*q_x);
           2*(q_x*q_z-q_w*q_y), 2*(q_y*q_z+q_w*q_x), q_w^2-q_x^2-q_y^2+q_z^2];


end