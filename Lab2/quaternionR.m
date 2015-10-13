function [R, q] = quaternionR(n, alpha)

q = [cos(alpha/2); n*sin(alpha/2)];
norm(q);

s = q(1);
v = q(2:4);

R = v*v' + s^2*eye(3) + 2*s*liu_crossop(v) + liu_crossop(v)^2;

end

