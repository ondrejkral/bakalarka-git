function v = ilspencskalna( A, b )
%ILSPENCSKALNA Proposition 3.56 in Hladik (?).

I = eye(dim(A));
R = mag(I - A);

v = infsup(-1,1)*(verifylss((I - R),mag(b)));

end

