function v = verifyinvert( A, b, p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

radiusVector = ilspencradius(p);
Acenter = ilspencmatrixcenter(A,p);
Acenterinv = inv(Acenter);
dimensions = ilspencmatrixdim(A);
M = intval(zeros(dimensions));
for k = 1:length(p)
   M = M + radiusVector(k)*abs(Acenterinv*intval(ilspencgetak(A,k)));
end

if verspectrad(sup(M)) < 1 
    v = 1;
else 
    v = 0;
end
end

