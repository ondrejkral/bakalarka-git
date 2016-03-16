function n2=vernorm2(A)
%    VERNORM2        Verified 2-norm of a complex (or real) matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a rectangular complex (or real) matrix A,
%        n2=vernorm2(A)
%    computes a verified enclosure n2 of the 2-norm of A (i.e., of norm(A,2)),
%    or fails (yields an interval of NaN's).
%    
%    See also VERSINGVAL, NORM.

%    Copyright 2008 Jiri Rohn.
%
%    Computed as the largest singular value of A using VERSINGVAL. 
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
n2=infsup(NaN,NaN);
if ~(nargin==1&&nargout<=1&&isreal(A)&&~isintval(A)) % wrong data
    setround(gr); return
end
sigma=versingval(A); % main part % interval vector of singular values
if isnan(sigma.inf(1)) % singular values not computed
    setround(gr); return 
end
% singular values computed
n2=sigma(1);
setround(gr);
