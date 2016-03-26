function [reg,As]=verregsing(A,t)
%    VERREGSING     Verified regularity/singularity of an interval matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a square interval matrix A,
%    [reg,As]=verregsing(A)
%    verifies regularity or singularity of A, or yields no verified result:
%
%    reg= 1           A verified regular,
%    reg= 0           A verified singular; As is a very tight ("almost thin") 
%                     interval matrix which is a part of A and is verified
%                     to contain a singular real matrix,
%    reg=-1           no verified result (A may be not square or not real).
%
%    [reg,As]=verregsing(A,1) [i.e., with additional input argument "1"]
%                     is the same as before, but it also produces screen
%                     output of the form "Expected remaining time: ... sec."
%                     This feature, however, slows down the actual computation.
%
%    Computational experience shows that checking singularity (reg=0) is
%    usually very fast; however, checking regularity (reg=1) occasionally may 
%    last long since the problem is NP-hard.
%
%    See also  VERIFYLSS, VERINTERVALHULL.

%    Copyright 2007 Jiri Rohn
%
%    Based on the algorithm "regularity" described in
%    J. Rohn, A Handbook of Results on Interval Linear Problems,
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
% set rounding to nearest
gr=getround;
setround(0);
% setting default output
[m,n]=size(A);
reg=-1; As=repmat(infsup(NaN,NaN),m,n); 
% checking compatibility of data
if (m~=n)||~isreal(A)
    % error('matrix not square or not real'),
    setround(gr); return
end
% creating full matrix of type intval
if ~isintval(A)
    A=infsup(A,A);                                    % allows for real input
end 
if issparse(A)
    A=full(A);                                        % sparse matrices not implemented yet
end
% time display
if (nargin==2)&&isequal(t,1)                          % t==1: display remaining time
    time=1;
else
    time=0;
end
% verified midpoint and radius
[ac,Delta]=vermidrad(A);                              
aci=inv(ac);
I=eye(n,n); I=infsup(I,I);
e=ones(n,1); e=infsup(e,e);
% checking regularity
% first, via verifylss
x=verifylss(A,e);
if ~isnan(x.inf(1))                                   % solution set bounded
    reg=1;                                            % regular; routine verifylss due to Rump
    setround(gr); return                              
end
% second, via strong regularity
M=inv(I-abs(aci)*Delta);
if isa(M,'double')                                    % bridging the current gap in verifylss
    M=infsup(M,M); 
end
if (~isnan(M.inf(1,1)))&&(all(all(M.inf>=0)))         % M verified nonnegative
    reg=1;                                            % strongly regular, Beeck 1975
    setround(gr); return                                     
end
% third, via positive definiteness                              
if isspd(ac'*ac-norm(Delta'*Delta,1)*I)               % verified positive definite
   reg=1;                                             % regular; Rex and Rohn, SIMAX 1998
   setround(gr); return                                    
end
% fourth, general regularity check
% step A: finding an appropriate right-hand side
b=ones(n,1);
acim=mid(aci);
if ~isnan(acim(1,1))                                  % heuristic procedure for finding appropriate b
    x=acim*b;                                         % performed in floating point arithmetic
    g=min(abs(x));                                    % Jansson and Rohn, SIMAX 1999
    for i=1:n
        for j=1:n
            xp=x-2*b(j)*acim(:,j);                    % corresponds to b(j)=-b(j);
            if min(abs(xp))>g,
                g=min(abs(xp)); x=xp; b(j)=-b(j);     % min(abs(xp)) increased
            end
        end
    end
end
b=infsup(b,b);
% step B: the regularity check itself
[x,AAs]=intervalhull(A,b,time);                       % main part: general check
if ~isempty(x)
    reg=1;                                            % verified regular
    setround(gr); return
end
if ~isempty(AAs)
    reg=0; As=AAs;                                    % verified singular
    setround(gr); return
end
reg=-1; As=repmat(infsup(NaN,NaN),n,n); 
setround(gr);
%
function [Ac,Delta]=vermidrad(A)
% computes verified midpoint and radius of A
% Ac, Delta are intval quantities
if ~isintval(A)
    Ac=infsup(A,A);
    Delta=infsup(zeros(size(A)),zeros(size(A)));
else
    Al=infsup(A.inf,A.inf);
    Au=infsup(A.sup,A.sup);
    Ac=   (Al+Au)/2;        
    Delta=(Au-Al)/2;        
end
