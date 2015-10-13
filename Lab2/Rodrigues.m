function [R] = Rodrigues( n,angle )
%RODRIGUEZ Calculates the rotation matrix

[p,q]=createONbasis(n);

nx = q*p'-p*q';

R = eye(3)+(1-cos(angle))*nx*nx+sin(angle)*nx;

end

