function [L,X]=verspectdec(A)
%    VERSPECTDEC        Verified spectral decomposition of a symmetric real matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For an n-by-n symmetric real matrix A,
%        [L,X]=verspectdec(A)
%    computes a diagonal n-by-n interval matrix L and an n-by-n
%    interval matrix X verified to contain a diagonal matrix Lo and an
%    orthogonal matrix Xo such that diag(Lo) is the vector of eigenvalues
%    of A and
%        A=Xo*Lo*Xo'
%    holds (a spectral decomposition of A). The entries of the interval
%    vector diag(L) are ordered in decreasing order. If no verified result
%    is achieved, then L, X consist of NaN's (L diagonal).
%    
%    See also VEREIGSYM, EIG.

%    Copyright 2008 Jiri Rohn.
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
%    2008-02-11   first version
%    2008-02-21   version for posting (p-coded)
%    2008-04-01   version for posting (decoded)
%
gr=getround;
setround(0);
[m,n]=size(A);
L=diag(repmat(infsup(NaN,NaN),n,1)); X=repmat(infsup(NaN,NaN),n,n); % setting default output
if (m~=n)||~isreal(A)||~isequal(A,A')||isintval(A) % wrong data
    setround(gr); return
end
if issparse(A)
    A=full(A); % creating full matrix; sparse matrices not implemented yet
end
A1=jxj(A); % recasting 
[L1,X1]=ol(A1); % main part
if isnan(L1.inf(1,1)) % no enclosures
    setround(gr); return
end
L1=jxj(L1); % recasting; ordered from largest to smallest
X1=jxj(X1); % recasting
if alldisjoint(diag(L1))==0 % eigenvalue enclosures not mutually disjoint
    setround(gr); return
end
for j=1:n % normalizing
    alpha=norm(X1(:,j),2);
    if alpha.inf==0 % normalizing not possible
        setround(gr); return
    end
    X1(:,j)=X1(:,j)/alpha;
end
L=L1;
X=X1;
setround(gr);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions jxj, alldisjoint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Y=jxj(X)
%Y=J*X*J
mn=size(X,1);
Y=X(:,mn:-1:1);
Y=Y(mn:-1:1,:);
%
function ad=alldisjoint(a)
% ad=1 if all entries of a disjoint, ad=0 otherwise
n=length(a);
for i=1:n-1
    for j=i+1:n
        if ~isnan(intersect(a(i),a(j)))
            ad=0; return
        end
    end
end
ad=1;

