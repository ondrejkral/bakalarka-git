function v = ilspencbauerskeel( A, b, p )
%ILSPENC Bauer-Skeel bounds to parametric interval systems.
%   Parametric generalization by Hladik (2012).

% Inicialization of general variables.
dimensions = ilspencmatrixdim(A);
M = intval(zeros(dimensions));
C = intval(zeros(dimensions(1),1));
I = eye(dimensions);

radiusVector = ilspencradius(p);
% precondition matrix
Acenter = ilspencmatrixcenter(A,p);
Acenterinv = inv(Acenter);
% aproximate solution
x = Acenterinv*ilspencbcenter(b ,p);

% Matrix M from Theorem 4.
for k = 1:length(p)
   M = M + radiusVector(k)*abs(Acenterinv*intval(ilspencgetak(A,k)));
end

% Summation in interval enclosure from Theorem 4.
for k = 1:length(p)
    C = C + radiusVector(k)*abs(Acenterinv*(ilspencgetak(A,k)*intval(x) - ilspencgetbk(b,k)));
end

s = verifylss(I - M,C);

v = hull(x - s, x + s);
end

