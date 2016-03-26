function v = verifyinvert2( A, b, p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[Ares, bres] = ilspencresidualform(A, b, p);

R = inv(mid(Ares));
C = eye(dim(Ares))-R*intval(Ares);

if verspectrad(sup(abs(C))) < 1 
    v = 1;
else 
    v = 0;

end

