function [xprim] = sandwichproduct(q, x0)

s = q(1);
v = q(2:4);

xprim = (v*v' + s^2*eye(3) + 2*s*liu_crossop(v) + liu_crossop(v)^2)*x0;

end

