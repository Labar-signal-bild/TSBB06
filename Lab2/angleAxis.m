function [ n, angle ] = angleAxis( R )
%ANGLEAXIS calculates the angle and axis from a rotation matirx

v = 0.5.*[R(3,2)-R(2,3); %Rx
          R(1,3)-R(3,1);
          R(2,1)-R(1,2)];
       
c = (trace(R)-1)/2;
s = norm(v);

if s = 0
    if c > 0
        n = s^2;
        angle = 0;
    else
        n = v/s;
        angle = pi;
    end
else
    n = v/s;
    angle =atan(s/c);
end
end

