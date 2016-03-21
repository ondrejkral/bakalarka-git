function v = verifyinvert3( A, b, p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[Ares, bres] = ilspencresidualform(A, b, p);

I = eye(dim(Ares));
R = mag(I - Ares);

if verspectrad(sup(R)) < 1 
    v = 1;
else 
    v = 0;

end
