function v = ilspencrump( A, b, iterations)
%ILSPENCRUMP Bounds to parametric interval systems by Rump(2010).

v = NaN;
% approximate inverse
R = inv(mid(A)); 
% % approximate solution
x = R*mid(b);
% iteration matrix
C = eye(dim(A))-R*intval(A);
Z = R*(b-A*intval(x));
X = Z; iter = 0;
while iter < iterations
    iter = iter+1;
    Y = X*infsup(0.9,1.1) + 1e-20*infsup(-1,1); % epsilon-inflation
    X = Z+C*Y; % interval iteration

    if all(in0(X,Y)), v = x + X; 
        return;
    end
end
end

