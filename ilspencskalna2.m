function v = ilspencskalna2( A, b, p )
%ILSPENCSKALNA2 Enclousure by Skalna (2006/4/1).
%   Enclosure of parametric linear equations by Iwona Skalna.

% Inicialization of general variables.
dimensions = ilspencmatrixdim(A);
Z = intval(zeros(dimensions,1));
D = intval(zeros(dimensions));

% precondition
Acenterinv = inv(ilspencmatrixcenter(A,p));

% aproximate solution
x = Acenterinv*ilspencbcenter(b ,p);

for k = 1:length(p)
   Z = Z + (ilspencgetbk(b,k)- ilspencgetak(A,k)*x)*p(k);
end
Z = Acenterinv*Z;

for k = 1:length(p)
   D = D + ilspencgetak(A,k)*p(k);
end
D = Acenterinv*D;

v = x + infsup(-1,1)*inv(mig(D))*mag(Z);
end

