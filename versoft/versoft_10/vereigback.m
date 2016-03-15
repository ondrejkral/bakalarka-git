function [lambda,X,ep]=vereigback(A)
%    VEREIGBACK        Verified backward error analysis of eigenpairs.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a square complex (or real) matrix A,
%    [lambda,X,ep]=vereigback(A)
%    computes a vector of eigenvalues lambda and a matrix of eigenvectors X
%    in the usual MATLAB way
%        [X,L]=eig(A); lambda=diag(L);
%    and additionally a vector ep with the following property: for each i
%    there exists a matrix, say A[i], verified to satisfy
%        max(max(abs(A-A[i]))) <= ep(i)
%    such that (lambda(i), X(:,i)) is verified to be an EXACT eigenpair
%    of A[i]. If A and lambda(i), X(:,i) are real then A[i] can be taken
%    real, otherwise it is complex in general. The maximal value of ep(i)
%    is usually very small (of order 1e-013 to 1e-016), which shows that
%    MATLAB computes eigenvalues and eigenvectors with great accuracy.
%
%    EXAMPLES. For a randomly generated complex 500x500 matrix we obtain 
%    >> n=500; rand('state',1); A=2*rand(n,n)-1+i*(2*rand(n,n)-1); [lambda,X,ep]=vereigback(A); epmax=max(ep)
%    epmax =
%      2.5875e-013.
%    On the other hand, for a "bad" 5x5 matrix gallery(5) we get
%    >> A=gallery(5); [lambda,X,ep]=vereigback(A); epmax=max(ep)
%    epmax =
%      2.7617e-011.

%    Copyright 2008 Jiri Rohn.
%
%    Based on the inequality (3.13) in
%    J. Rohn, A Handbook of Results on Interval Linear Problems,
%    posted at http://www.cs.cas.cz/~rohn,
%    which also holds for complex eigenpairs (unpublished).
%
%    This work was supported by the Czech Republic National Research
%    Program "Information Society", project 1ET400300415. 
%
%    WARRANTY
%
%    Because the program is licensed free of charge, there is 
%    no warranty for the program, to the extent permitted by applicable
%    law. Except when otherwise stated in writing the copyright holder
%    and/or other parties provide the program "as is" without warranty
%    of any kind, either expressed or implied, including, but not
%    limited to, the implied warranties of merchantability and fitness
%    for a particular purpose. The entire risk as to the quality and
%    performance of the program is with you. Should the program prove
%    defective, you assume the cost of all necessary servicing, repair
%    or correction.
%
gr=getround;
setround(0);
[m,n]=size(A);
lambda=repmat(NaN,n,1); X=repmat(NaN,m,n); ep=lambda;
if m~=n
    setround(gr); return
end
I=eye(n,n); 
[X,L]=eig(A); % A*X=X*L
lambda=diag(L);
for i=1:n
    ll=lambda(i);
    ll=infsup(ll,ll);
    xx=X(:,i);
    xx=infsup(xx,xx);
    epi=norm((A-ll*I)*xx,'inf')/norm(xx,1); % main formula
    ep(i)=epi.sup;
end
setround(gr);
