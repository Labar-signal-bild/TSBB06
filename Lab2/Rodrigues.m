function [R] = Rodrigues( n,angle )
%RODRIGUEZ Calculates the rotation matrix

p = [-n(3);0;n(1)];
q = cross(p,n);

p = p/norm(p);
q = q/norm(q);

nx = q*p'-p*q';


R = eye(3)+(1-cos(angle)).*nx^2+sin(angle).*nx;

end

