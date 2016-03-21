function lam=vereigsym(A)
%    VEREIGSYM     Verified eigenvalues of a symmetric interval (or real) matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    CONVENTION. A real symmetric n-by-n matrix Ao has n real eigenvalues
%    lambda(Ao,i), i=1,...,n. We order them in nonincreasing sequence, i.e., to satisfy
%        lambda(Ao,1) >= lambda(Ao,2) >= ... >= lambda(Ao,n).
%        
%    For a symmetric interval matrix A (that is, satisfying A'=A),
%        lam=vereigsym(A)
%    produces an interval vector lam verified to satisfy
%        lambda(Ao,i) in lam(i), i=1,...,n
%    for each SYMMETRIC Ao in A. Multiplicity of eigenvalues is taken into
%    account. Nothing is said of eigenvalues of nonsymmetric matrices in A.
%    If no verified result if given, then lam is an interval vector of NaN's.
%
%    The interval vector lam has an additional property: for each i and j,
%    the intervals lam(i) and lam(j) are either identical, or disjoint.
% 
%    Accordingly, for a real symmetric matrix A,   
%        lam=vereigsym(A)
%    produces an interval vector lam verified to satisfy
%        lambda(A,i) in lam(i), i=1,...,n,
%    i.e., enclosing the vector of eigenvalues of A.

%    Copyright 2008 Jiri Rohn
%
%    Based partially on the section "Eigenvalues of symmetric matrices" in
%    J. Rohn, A Handbook of Results on Interval Linear Problems,
%    posted at http://www.cs.cas.cz/~rohn
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
%    2008-02-02   first version
%    2008-02-08   version for posting (p-coded)
%    2008-04-01   version for posting (decoded)
%
gr=getround;
setround(0);
% setting default output
[m,n]=size(A);
lam=repmat(infsup(NaN,NaN),n,1);
if (m~=n)||~isreal(A)||~isequal(A,A') % wrong data
    setround(gr); return
end
% creating full matrix 
if issparse(A)
    A=full(A); % sparse matrices not implemented yet
end
% case of real matrix
if ~isintval(A)
    L=ol(A); % diagonal matrix
    lam=diag(L); % vector of verified eigenvalues of A
    lam=lam(n:-1:1); % reorder from largest to smallest
    setround(gr); return
end
% case of interval matrix
Ac=mid(A); Delta=rad(A); 
if ~(isequal(Ac,Ac')&&isequal(Delta,Delta')&&all(all(in(A,midrad(Ac,Delta)))))
    setround(gr); return
end
L=ol(Ac); % diagonal matrix
if isnan(L.inf(1,1)) % verified eigenvalues of Ac not computed
    setround(gr); return
end
lambda=diag(L); % vector of verified eigenvalues of Ac
rho=verspectrad(Delta); % verified spectral radius
if isnan(rho.inf) % if not computed, replaced by norm
    rho=norm(infsup(Delta,Delta),1);
end
lam=lambda+midrad(zeros(n,1),rho.sup*ones(n,1)); % main formula % Handbook, p. 29, (3.14)
lam=lam(n:-1:1); % reorder from largest to smallest
setround(gr);
