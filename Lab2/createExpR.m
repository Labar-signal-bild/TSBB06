function [R] = createExpR(n, alpha)

[p,q] = createONbasis(n);
E = [n (p+i*q)/sqrt(2) (p-i*q)/sqrt(2)];
D = [1 0 0; 0 exp(i*alpha) 0; 0 0 exp(-i*alpha)];

R = E*D*conj(E)';

end

