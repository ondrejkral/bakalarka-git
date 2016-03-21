function dt=verdet(A)
%    VERDET        Verified determinant of a complex (or real) matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a square complex (or real) matrix A,
%        dt=verdet(A)
%    computes a verified enclosure dt of the determinant of A,
%    or fails (yields an interval of NaN's).
%
%    For an integer matrix A of size up to 7, the determinant value may
%    turn out to be exact.
%
%    EXAMPLE. 
%    A =
%       -56  -100    74    48   -69   -29   -88
%        45   -33   -95   -31   -17    49     8
%       -86   -45    45    77   -81    30    -9
%        93   -91    70   -31   -10    88    73
%       -58   -81    46   -88    74    67    71
%       -68   -18    91    44   -22    -6    -6
%        28    63    31    92   -49    26    57
%    >> format long, dt=verdet(A), width=dt.sup-dt.inf
%    intval dt = 
%      1.0e+012 *
%    [   2.03813924493800,   2.03813924493800] 
%    width =
%         0
%    As we can see, the resulting value of the determinant is huge, but
%    still computed exactly.
%
%    See also DET.

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
%    2008-02-05   first version
%    2008-02-14   version for posting
%    2008-03-11   case of integer matrix added
%    2008-03-31   version for posting
%
gr=getround;
setround(0);
[m,n]=size(A);
dt=infsup(NaN,NaN);
if ~(nargin==1&&nargout<=1&&m==n&&~isintval(A)) % wrong data
    setround(gr); return
end
L=ol(A); 
l=diag(L); 
if isnan(l.inf(1)) % not computed
    setround(gr); return 
end
% computed
dt=prod(l);
if isreal(A)
    dt=real(dt);
    % case of integer matrix
    if isequal(round(A),A)
        dt=infsup(ceil(dt.inf),floor(dt.sup));
    end
end
setround(gr);
