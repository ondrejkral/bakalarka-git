function [hs,As]=verhurwstab(A,t)
%    VERHURWSTAB    Verified Hurwitz stability of an interval matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a square interval matrix A,
%    [hs,As]=verhurwstab(A)
%    verifies Hurwitz stability or not-Hurwitz-stability of A,
%    or yields no verified result:
%
%    hs= 1           A verified Hurwitz stable,
%    hs= 0           A verified not to be Hurwitz stable; As is a very
%                    tight ("almost thin") interval matrix which is a part
%                    of A and is verified to contain a not-Hurwitz-stable
%                    real matrix,    
%    hs=-1           no verified result.
%
%    [hs,As]=verhurwstab(A,1) [i.e., with additional input argument "1"]
%                    is the same as before, but it also produces screen
%                    output of the form "Expected remaining time: ... sec."
%                    This feature, however, slows down the actual computation.
%
%    The algorithm employs only a sufficient (not necessary) condition and,
%    unfortunately, the output hs=-1 is rather frequent. The case hs=0
%    occurs only for symmetric interval matrices.
%
%    See also ISSPD, VERPOSDEF.

%    Copyright 2007 Jiri Rohn
%
%    Based on the algorithm "hurwitzstab" described in
%    J. Rohn, A handbook of results on interval linear problems,
%    posted at http://www.cs.cas.cz/~rohn
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
% setting default output
[m,n]=size(A);
hs=-1; As=repmat(infsup(NaN,NaN),m,n);
% checking compatibility of data
if (m~=n)||~isreal(A)                            % error('matrix not square or not real')
    setround(gr); return
end
% creating full matrix of type intval
if ~isintval(A)
    A=infsup(A,A);                               % allows for real input
end 
if issparse(A)
    A=full(A);                                   % sparse matrices not implemented yet
end
% time display
if (nargin==2)&&isequal(t,1)                     % t==1: display remaining time          
    time=1;
else
    time=0;
end
% symmetrization
if isequal(A,A')                                 
    Asym=A;
else
    Asym=(A+A')/2;                               
end
% main check
[pd,AAs]=verposdef(-Asym,time);                  % positive definiteness check of -Asym; Rohn, SIMAX 1994, Thm. 6
% output cases
switch pd
    case  1 
        hs=1;                                    % verified Hurwitz stable; Rohn, SIMAX 1994, Thm. 7
    case  0
        if isequal(A,A')
            hs=0; As=-AAs;                       % verified not Hurwitz stable; As of type intval
        else
            hs=-1;                               % no verified result
        end
    case -1
        hs=-1;                                   % no verified result
end
setround(gr);
