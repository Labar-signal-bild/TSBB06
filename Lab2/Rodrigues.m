function [R] = Rodrigues( n,angle )
%RODRIGUEZ Calculates the rotation matrix

q = [0;-n(3);n(2)];
p = cross(n,q);
n = n'./norm(n);
q = q./norm(q);
p = p'./norm(p);

nx = q*p'-p*q'; 

R = eye(3)+(1-cos(angle)).*nx^2+sin(angle).*nx;

end

