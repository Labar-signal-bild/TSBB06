function [p, q] = createONbasis(n)

p = [-n(3);0;n(1)];
q = cross(p,n);

p = p/norm(p);
q = q/norm(q);

end
