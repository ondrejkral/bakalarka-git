function v = ilspenchbrref( A,b,p,x)
%ILSPENCHBRREF Refinement of Hans-Bliek-Rohn generalization.
% Refinement by Hladik (2012).

% p = vector of interval values in which parameters lies
% x = enclosure obtained from Hans-Bliek-Rohn method

% Inicialization of general variables.
dimensions = ilspencmatrixdim(A);
Y = intval(zeros(dimensions)); y = intval(zeros(dimensions(1),1));
Z = intval(zeros(dimensions)); z = intval(zeros(dimensions(1),1));
I = eye(dimensions);
v = intval(zeros(dimensions(1),1));

radiusvector = ilspencradius(p);
% Precondition matrix.
Acenter = ilspencmatrixcenter(A,p);
Acenterinv = inv(Acenter);
% approximate solution 
x1 = Acenterinv*ilspencbcenter(b ,p);

for k = 1:length(p)
    a = Acenterinv*(ilspencgetak(A,k)*x - ilspencgetbk(b,k));
    for j = 1:dimensions(1)
        if inf(a(j)) >= 0
            Y(j,:) = Y(j,:) + radiusvector(k)*Acenterinv(j,:)*intval(ilspencgetak(A,k));
            y(j) = y(j) + radiusvector(k)*Acenterinv(j,:)*intval(ilspencgetbk(b,k));
        elseif sup(a(j)) <= 0
            Y(j,:) = Y(j,:) - radiusvector(k)*Acenterinv(j,:)*intval(ilspencgetak(A,k));
            y(j) = y(j) - radiusvector(k)*Acenterinv(j,:)*intval(ilspencgetbk(b,k));
        else
            Z(j,:) = Z(j,:) + radiusvector(k)*abs(Acenterinv(j,:)*intval(ilspencgetak(A,k)));
            z(j) = z(j) + radiusvector(k)*abs(Acenterinv(j,:)*intval(ilspencgetbk(b,k)));
        end
    end
end

% M-asteriks from refinement algorithm 2.
M0 = I - abs(Y) - Z;
M = inv(M0);
% x upper-index zero from refinement algorithm 2.
x0 = verifylss(M0,abs(x1) - y + z);

for i = 1:dimensions(1)
    u = x0(i) + (x1(i) - abs(x1(i)))*M(i,i);
    d = - x0(i) + (x1(i) + abs(x1(i)))*M(i,i);
    coef = 1/(2*M(i,i) - 1);
    v(i) = hull(min(inf(d),inf(coef*d)),max(sup(u),sup(coef*u)));
end
end

