function [q] = rotvec2quat(rotvec)
%Convert a rotation vector to a quaternion
%   Input: 
%           rotvec: 1x3 vector
%   Output: 
%                q: 1x4 vector     
    assert(all(size(rotvec) == [1, 3]));

    phi = norm(rotvec);
    if phi > 1e-8
        u = rotvec / phi;
        q = [cos(phi/2) u * sin(phi / 2)];
    else
        q = [1 0 0 0];
    end

end

