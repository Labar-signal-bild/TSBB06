function [ n, angle ] = angleAxis( R )
%ANGLEAXIS calculates the angle and axis from a rotation matirx

v = 0.5.*[R(3,2)-R(2,3); %Rx
          R(1,3)-R(3,1);
          R(2,1)-R(1,2)];
       
c = (trace(R)-1)/2;
s = norm(v);

if s == 0
    display('---- s=0 ----');
    if c > 0
        display('---- c>0 ----');
        n = 2*rand(3,1)-1;
        n = n/norm(n);
        angle = 0;
    else
        display('---- c<=0 ----');
        [evec eval] = eig(R+eye(3));
        pos = find(max(eval)==eval);
        n = evec(:,pos);
        
        %n = v/s;
        angle = pi;
        
    end
else
    %display('----normal case ----')
    n = v/s;
    angle = atan(s/c);
end
end

