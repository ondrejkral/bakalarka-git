function v = ilspencbsref( A,b,p,x)
%ILSPENCBSREF Refinement of Bauer-Skeel generalization.
%   Refinement by Hladik (2012)

% p = vector of interval values in which parameters lies
% x = enclosure obtained from Bauer-Skeel method

% Inicialization of general variables.
dimensions = ilspencmatrixdim(A);
Y = intval(zeros(dimensions)); y = intval(zeros(dimensions(1),1));
Z = intval(zeros(dimensions)); z = intval(zeros(dimensions(1),1));
I = eye(dimensions);

radiusvector = ilspencradius(p);
% Precondition matrix.
Acenter = ilspencmatrixcenter(A,p);
Acenterinv = inv(Acenter);
% x-asterisk from Theorem 4.
x1 = intval(Acenterinv*ilspencbcenter(b ,p));

for k = 1:length(p)
    a = Acenterinv*(ilspencgetak(A,k)*x - ilspencgetbk(b,k));
    for j = 1:dimensions(1)
        if inf(a(j)) >= 0
            Y(j,:) = Y(j,:) + radiusvector(k)*Acenterinv(j,:)*intval(ilspencgetak(A,k));
            y(j) = y(j) + radiusvector(k)*Acenterinv(j,:)*(ilspencgetak(A,k)*x1 - ilspencgetbk(b,k));
        elseif sup(a(j)) <= 0
            Y(j,:) = Y(j,:) - radiusvector(k)*Acenterinv(j,:)*intval(ilspencgetak(A,k));
            y(j) = y(j) - radiusvector(k)*Acenterinv(j,:)*(ilspencgetak(A,k)*x1 - ilspencgetbk(b,k));
        else
            Z(j,:) = Z(j,:) + radiusvector(k)*abs(Acenterinv(j,:)*intval(ilspencgetak(A,k)));
            z(j) = z(j) + radiusvector(k)*abs(Acenterinv(j,:)*(ilspencgetak(A,k)*x1 - ilspencgetbk(b,k)));
        end
    end
end

refinement = verifylss(I - abs(Y) - Z, y + z);

v = hull(x1 - refinement, x1 + refinement);
end

