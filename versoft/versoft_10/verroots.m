function zs=verroots(a)
%    VERROOTS       Verified roots of a complex (or real) polynomial.   
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a complex (or real) vector a of length n+1 with a(1)=1, 
%        zs=verroots(a)
%    returns a vector zs of verified ALL roots of the polynomial
%        a(1)*z^n+...+a(n)*z+a(n+1),   (1)
%    or yields no verified result.
%
%    COMMENT. It must be a(1)=1 (otherwise the computation breaks down), so
%    that the polynomial is monic.
%
%    Possible outcomes:
%
%    ~isnan(zs.inf(1,1)) :  zs is  a (generally complex) interval vector
%                           verified to contain the vector of all roots of
%                           the polynomial (1); it is sorted in ascending
%                           order of the real parts of the midpoints,   
%     isnan(zs.inf(1,1)) :  no verified result (the interval vector zs
%                           consists of NaN's).  
%
%    See also VEREIG, VERIFYEIG, EIG, ROOTS.

%    Copyright 2008 Jiri Rohn
%
%    Polynomial roots computed as eigenvalues of the companion matrix.
%
%    This work was done during author's employment at the Anglo-American
%    University in Prague, Czech Republic.
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
%    2008-05-15   first version 
%    2008-05-28   version for posting
%
gr=getround;
setround(0);
% a: coefficients of the monic polynomial a(1)*z^n+...+a(n)*z+a(n+1), a(1)=1
a=a(:);
n=length(a)-1;  % degree of the polynomial
if n<1||a(1)~=1 % must be monic
    zs=repmat(infsup(NaN,NaN),n,1);
    setround(gr); return
end
% companion matrix
C=diag(ones(1,n-1),-1); 
C(1,:)=-a(2:n+1); % first row
% verified eigenvalues of C
L=ol(C);    % diagonal matrix 
zs=diag(L); % vector of verified zeros
% sorting in ascending order of real parts of the midpoints
[y,I]=sort(real(mid(zs))); 
zs=zs(I);
setround(gr);