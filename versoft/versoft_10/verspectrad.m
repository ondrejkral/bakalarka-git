function rho=verspectrad(A)
%    VERSPECTRAD        Verified spectral radius of a complex (or real) matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a square complex (or real) matrix A,
%        rho=verspectrad(A)
%    computes a verified enclosure rho of the spectral radius of A,
%    or fails (yields an interval of NaN's).
%    
%    See also EIG.

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
%    2008-02-01   first version
%    2008-02-21   version for posting (p-coded)
%    2008-04-01   version for posting (decoded)
%
gr=getround;
setround(0);
[m,n]=size(A);
rho=infsup(NaN,NaN);
if ~(nargin==1&&nargout<=1&&m==n&&~isintval(A)) % wrong data
    setround(gr); return
end
L=ol(A); % main part
lambda=diag(L); 
if isnan(lambda.inf(1)) % eigenvalues not computed
    setround(gr); return 
end
% eigenvalues computed
a=abs(lambda);
[rhosup,i]=max(a.sup); % upper bound
rhoinf=a.inf(i); % current lower bound
for j=1:n
    if ~isnan(intersect(a(i),a(j)))
        b=a(j);
        rhoinf=min(b.inf,rhoinf); % update of the lower bound
    end
end
rho=infsup(rhoinf,rhosup); % verified spectral radius
setround(gr);
