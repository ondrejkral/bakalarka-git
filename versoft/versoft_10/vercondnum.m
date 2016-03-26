function cn=vercondnum(A)
%    VERCONDNUM        Verified 2-norm condition number of a square real matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a square real matrix A,
%        cn=vercondnum(A)
%    computes a verified 2-norm condition number cn of a square real matrix A
%    (i.e., norm(A,2)*norm(inv(A),2)), or fails (yields an interval of NaN's).  
%
%    See also VERSINGVAL, COND.

%    Copyright 2008 Jiri Rohn.
%
%    Computed as the ratio of the verified largest and smallest singular
%    value of A using VERSINGVAL.
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
%    2008-02-05   first version
%    2008-02-14   version for posting
%
gr=getround;
setround(0);
[m,n]=size(A);
cn=infsup(NaN,NaN);
if ~(nargin==1&&nargout<=1&&m==n&&isreal(A)&&~isintval(A)) % wrong data
    setround(gr); return
end
sigma=versingval(A); % main part % interval vector of singular values
if isnan(sigma.inf(1)) % singular values not computed
    setround(gr); return 
end
% singular values computed
if isequal(sigma.inf(n),0) % lower bound on sigma(n) is zero
    setround(gr); return
end
cn=sigma(1)/sigma(n); % main formula
setround(gr);
