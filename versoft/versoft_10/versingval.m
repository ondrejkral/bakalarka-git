function sig=versingval(A)
%    VERSINGVAL     Verified singular values of a rectangular interval (or real) matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    REMINDER. A real m-by-n matrix Ao has q=min(m,n) singular values
%       sigma(Ao,1) >= sigma(Ao,2) >= ... >= sigma(Ao,q) >= 0.
%        
%    For a rectangular m-by-n interval matrix A 
%        sig=versingval(A)
%    produces a nonnegative interval vector sig of length q=min(m,n) verified to satisfy
%        sigma(Ao,i) in sig(i), i=1,...,q
%    for each Ao in A. If no verified result if given, then sig is an
%    interval vector of NaN's. 
%
%    Accordingly, for a real m-by-n matrix A,   
%        sig=versingval(A)
%    produces an interval vector sig verified to satisfy
%        sigma(A,i) in sig(i), i=1,...,q,
%    i.e., enclosing the vector of singular values of A.
%
%    See also VEREIGSYM.

%    Copyright 2008 Jiri Rohn
%
%    Based on the Jordan-Wielandt theorem, with the verified eigenvalues of
%    the respective symmetric matrix being computed by VEREIGSYM.
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
%    2008-02-09   first version
%    2008-02-10   version for posting
%
gr=getround;
setround(0);
[m,n]=size(A);
q=min(m,n);
sig=repmat(infsup(NaN,NaN),q,1); % setting default output
if ~isreal(A)
    setround(gr); return
end
if ~isintval(A) % construction of the Jordan-Wielandt matrix
    Ajw=[zeros(m,m) A; ctranspose(A) zeros(n,n)];
else
    Ajw=[infsup(zeros(m,m),zeros(m,m)) A; ctranspose(A) infsup(zeros(n,n),zeros(n,n))];
end
sig=vereigsym(Ajw); % eigenvalues of Ajw (m+n of them) in nonincreasing order
if isnan(sig.inf(1)) % not computed
    setround(gr); return
end
sig=sig(1:q); % reduction to first q 
o=zeros(q,1);
sig=infsup(max(sig.inf,o),max(sig.sup,o)); % ensuring nonnegativity of singular values
setround(gr);

