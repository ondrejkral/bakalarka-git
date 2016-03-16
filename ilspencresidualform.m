function [ Ares, bres] = ilspencresidualform( A, b, p )
%ILSPENCRESIDUALFORM Convert interval linear system to residual form.

dimensions = ilspencmatrixdim(A);
Ares = intval(zeros(dimensions));
bres = intval(zeros(dimensions(1),1));
% Precondition matrix.
C = ilspencmatrixcenter(A,p);
Cinv = inv(C);
% x-asterisk from Theorem 4.
x = verifylss(C,ilspencbcenter(b ,p));

for k = 1:length(p)
    Ares = Ares + p(k)*(Cinv*ilspencgetak(A,k));
end

for k = 1:length(p)
    bres = bres + p(k)*(Cinv*(ilspencgetbk(b,k) - ilspencgetak(A,k)*x));
end

end

