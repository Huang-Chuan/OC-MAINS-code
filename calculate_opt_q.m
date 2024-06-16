function [q_ast] = calculate_opt_q(q_km1, q_k, q_hat, g)
% calculate the optimal q_ast,  s.t. 1/2 * |q_ast - q_hat|_2^2
% subject to R{q_k} g = R{q_ast} R{q_km1} g
%            q_ast * q_ast' = 1

u1 = q2r(q_km1) * g;
u2 = q2r(q_k) * g;

xi_1= [0;u1+u2];
xi_2= [(u1+u2)'*(u1+u2); vect2skew(u1-u2)*(u1+u2)];

xi_1_normalized = xi_1 / norm(xi_1);
xi_2_normalized = xi_2 / norm(xi_2);

x_ast = q_hat*xi_1_normalized / sqrt((q_hat*xi_1_normalized)^2 + ...
    q_hat*xi_2_normalized)^2;
y_ast = q_hat*xi_2_normalized / sqrt((q_hat*xi_1_normalized)^2 + ...
    q_hat*xi_2_normalized)^2;

q_ast = x_ast * xi_1_normalized + y_ast * xi_2_normalized;


end