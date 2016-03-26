function [F,E]=vermatfun(f,A,r)
%    VERMATFUN        Verified matrix function of a complex (or real) matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a complex (or real) function f of one variable and for a square
%    complex (or real) matrix A, 
%        [F,E]=vermatfun(f,A)
%    computes a verified enclosure F of the matrix function f(A), or fails
%    (yields a matrix of NaN's). 
%
%    f must be given as an inline function (i.e., an expression between
%    apostrophes), as e.g. in
%        F=vermatfun('(1+x)/x',A).
%
%    The output F is typically a complex interval matrix even if the exact
%    result is real. If you know beforehand that the result is real (as it
%    is e.g. with the functions exp(x), sin(x), cos(x) etc. for real
%    argument), use
%        F=vermatfun(f,A,1)
%    for real output. 
%
%    The structured array E explains reasons for NaN output. It has three
%    fields: E.error, E.where, E.value. 
%
%    EXAMPLE (Golub and van Loan [1], p. 567). 
%    >> A =
%       -49    24
%       -64    31
%    >> format long, F=vermatfun('exp(x)',A,1)
%    intval F = 
%    [  -0.73575875814481,  -0.73575875814470] [   0.55181909965806,   0.55181909965814] 
%    [  -1.47151759908836,  -1.47151759908816] [   1.10363824071549,   1.10363824071565] 
%
%    [1] G. H. Golub and C. V. van Loan, Matrix Computations, The Johns
%    Hopkins University Press, Baltimore 1996.

%    Copyright 2008 Jiri Rohn.
%
%    Based of the fact that if 
%        A=X*L*inv(X) 
%    is an eigenvalue decomposition of a diagonalizable matrix A, then
%        f(A)=X*f(L)*inv(X)
%    (see [1], Corollary 11.1.2).
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
%    History
%
%    2008-03-23   first version
%    2008-03-31   output variable E added; version for posting
%
gr=getround;
setround(0);
[m,n]=size(A);
F=repmat(infsup(NaN,NaN),n,n);
E.error='vermatfun: none';
E.where='NaN';
E.value='NaN';
if ~(2<=nargin&&nargin<=3&&nargout<=2&&m==n&&~isintval(A)) % wrong data
    E.error='vermatfun: wrong data';
    setround(gr); return
end
if issparse(A)
    A=full(A); % sparse matrices not implemented yet
end
[L,X]=ol(A); % main part
if isnan(L.inf(1,1)) % no result
    E.error='vermatfun: ol fails: eigenvalue decomposition not computed';
    setround(gr); return
end
Y=inv(X);
if isnan(Y.inf(1,1)) % inverse not found
    E.error='vermatfun: inv fails: inverse not computed';
    setround(gr); return
end
% X, L, Y computed
f=inline(f);
LL=L; % preallocation
for j=1:n
    LL(j,j)=f(L(j,j)); % LL=f(L)
    if isinf(LL(j,j))||isnan(LL(j,j))
        E.error='vermatfun: function value not finite';
        E.where=['j = ',int2str(j)];
        E.value=['LL(j,j) = [',num2str(LL(j,j).inf),', ',num2str(LL(j,j).sup),']'];
        setround(gr); return
    end
    X(:,j)=LL(j,j)*X(:,j); % X=X*LL
end
F=X*Y; % F=X*LL*Y=X*f(L)*inv(X) (basic formula)
if isequal(nargin,3)&&isequal(r,1) % case of real result
    F=real(F);
end
setround(gr);
