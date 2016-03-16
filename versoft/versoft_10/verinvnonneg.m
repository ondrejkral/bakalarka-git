function [nonneg,As]=verinvnonneg(A)
%    VERINVNONNEG        Verified nonnegative invertibility of an interval matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a square interval (or real) matrix A,
%    [nonneg,As]=verinvnonneg(A)
%    verifies inverse nonnegativity of A, or not-inverse-nonnegativity of A,
%    or yields no verified result: 
%
%    nonneg= 1           A verified inverse nonnegative,
%    nonneg= 0           A verified not to be inverse nonnegative; As is a
%                        matrix in A (always one of the two bounds) which
%                        is verified not to be inverse nonnegative, 
%    nonneg=-1           no verified result.

%    Copyright 2008 Jiri Rohn
%
%    Based on the result by Kuttler, Math. of Computation 1971; see also
%    J. Rohn, A Handbook of Results on Interval Linear Problems,
%    posted at http://www.cs.cas.cz/~rohn, Section 3.9.
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
nonneg=-1; As=repmat(infsup(NaN,NaN),m,n);
if m~=n,                                                   % error('matrix not square')
    setround(gr); return
end
if ~isintval(A)
    A=infsup(A,A);                                         % allows for real input
end
Al=inf(A); Al=infsup(Al,Al);
Au=sup(A); Au=infsup(Au,Au);
Bl=inv(Al); 
Bu=inv(Au);
if (~isnan(Bl.inf(1,1)))&&(all(all(Bl.inf>=0)))&&(~isnan(Bu.inf(1,1)))&&(all(all(Bu.inf>=0)))
    nonneg=1;                                              % verified inverse nonnegative; Kuttler, Math. of Comp. 1971
    setround(gr); return
end
if (~isnan(Bl.inf(1,1)))&&(any(any(Bl.sup<0)))
    nonneg=0; As=Al;                                       % A.inf verified not inverse nonnegative
    setround(gr); return 
end
if (~isnan(Bu.inf(1,1)))&&(any(any(Bu.sup<0)))
    nonneg=0; As=Au;                                       % A.sup verified not inverse nonnegative
    setround(gr); return 
end
setround(gr);
