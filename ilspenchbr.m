function v = ilspenchbr( A, b, p )
%ILSPENCHBR Hansen-Bliek-Rohn bounds to parametric interval systems.
%   Parametric generalization by Hladik (2012).

% Inicialization of general variables.
dimensions = ilspencmatrixdim(A);
M = intval(zeros(dimensions));
C = intval(zeros(dimensions(1),1));
I = eye(dimensions);
v = intval(zeros(dimensions(1),1));

radiusVector = ilspencradius(p);
% Precondition matrix.
Acenter = ilspencmatrixcenter(A,p);
Acenterinv = inv(Acenter);
% x-asterisk from Theorem 4.
x = intval(Acenterinv*ilspencbcenter(b ,p));

% Matrix M from Theorem 4.
for k = 1:length(p)
   M = M + radiusVector(k)*abs(Acenterinv*intval(ilspencgetak(A,k)));
end

%M-asteriks from Theorem 5.
M0 = I - M;
M2 = inv(M0);

% Summation in x upper-index zero definition from Theorem 5.
for k = 1:length(p)
    C = C + radiusVector(k)*abs(Acenterinv*intval(ilspencgetbk(b,k)));
end

% x upper-index zero from Theorem 5. !!!
x0 = verifylss(M0,abs(x) + C);


for i = 1:length(x)
    u = x0(i) + (x(i) - abs(x(i)))*M2(i,i);
    d = -x0(i) + (x(i) + abs(x(i)))*M2(i,i);
    coef = 1/(2*M2(i,i) - 1);
    
    v(i) =  hull(min(inf(d),inf(coef*d)),max(sup(u),sup(coef*u)));
end

end

