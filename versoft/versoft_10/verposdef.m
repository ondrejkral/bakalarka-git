function [pd,As]=verposdef(A,t)
%    VERPOSDEF     Verified positive definiteness of an interval matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a symmetric interval matrix A,
%    [pd,As]=verposdef(A)
%    verifies positive definiteness or not-positive-definiteness of A, or
%    yields an unverified result: 
%
%    pd= 1           A verified positive definite,
%    pd= 0           A verified not to be positive definite; As is a very
%                    tight ("almost thin") interval matrix which is a part
%                    of A and is verified to contain a not-positive-definite  
%                    real matrix, 
%    pd=-1           no verified result (A may be not square or not real
%                    or not symmetric).
%
%    [pd,As]=verposdef(A,1) [i.e., with additional input argument "1"]
%                    is the same as before, but it also produces screen
%                    output of the form "Expected remaining time: ... sec."
%                    This feature, however, slows down the actual computation.
%
%    Computational experience shows that checking not-positive-definiteness (pd=0) 
%    is usually very fast; however, checking positive definiteness (pd=1)
%    occasionally may last long since the problem is NP-hard.
%
%    See also  ISSPD, VERREGSING.

%    Copyright 2007 Jiri Rohn
%
%    Based on the algorithm "posdefness" described in
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
pd=-1; As=repmat(infsup(NaN,NaN),m,n);
% checking compatibility of data
if (m~=n)||~isreal(A)||~isequal(A,A')
    % error('matrix not square or not real or not symmetric'),
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
% checking via isspd
if isspd(A)                                           % verified positive definite; isspd due to Rump
    pd=1; 
    setround(gr); return                              
end
% inspecting endpoint matrices
if ~isspd(A.inf)&&~isspd(A.sup)                       % none of the endpoints found positive definite
    i=verpd(A.inf);
    if i==0
        pd=0; As=infsup(A.inf,A.inf);                 % A.inf verified not positive definite; As of type intval for the sake of unified output
        setround(gr); return
    end
    s=verpd(A.sup);
    if s==0
        pd=0; As=infsup(A.sup,A.sup);                 % A.sup verified not positive definite: As of type intval for the sake of unified output
        setround(gr); return
    end
    if (i<0)&&(s<0)
        pd=-1;                                        % no result for both bounds
        setround(gr); return
    end                                               
end                                                   % at least one endpoint matrix verified positive definite
% checking via regularity                             % at this point, pos. definiteness of A is equivalent to regularity (Rohn, SIMAX 1994, Thm. 3) 
[reg,AAs]=verregsing(A,time);                         % main part: regularity check
switch reg
    case  1
        pd=1;                                         % verified positive definite
        setround(gr); return
    case  0
        AAs=(AAs+AAs')/2;                             % AAs singular but generally not symmetric
        AAs=intersect(AAs,A);
        if any(any(isnan(AAs)))                       % AAs not a part of A
            pd=-1; As=repmat(infsup(NaN,NaN),n,n);    % enclosure of a not-positive-definite matrix in A not found
            setround(gr); return
        end
        pd=0; As=AAs;                                 % AAs part of A, verified not positive definite (singular)
        setround(gr); return
    case -1
        pd=-1; As=repmat(infsup(NaN,NaN),n,n);        % no verified result
        setround(gr); return
end
%
% Subfunction
%
function [pd,x]=verpd(A)
%    VERPD    Verified positive definiteness of a real matrix.
%
%    For a symmetric real matrix A,
%    [pd,x]=verpd(A)
%    verifies positive definiteness or not-positive-definiteness of A,
%    or yields no verified result:
%
%    pd= 1 :  A verified positive definite,
%    pd= 0 :  A verified not positive definite, and
%             x is verified to satisfy x'*A*x<=0,
%    pd=-1 :  no verified result,
%    pd=-2 :  A not real or not symmetric.
%
%    If pd~=0, then x is a vector of NaN's.
%
%    See also ISSPD.

%    Copyright 2007 Jiri Rohn
%    The routine ISSPD copyrighted by Siegfried M. Rump.
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
% checking data
n=size(A,1);
pd=-1; x=repmat(NaN,n,1);
if isintval(A)||~isequal(A,A')
    pd=-2; 
    setround(gr); return
end
% verifying positive definiteness
if isspd(A)==1                                        % routine isspd by S. M. Rump
    pd=1;                                             % verified positive definite
    setround(gr); return
end
% verifying not-positive-definiteness
% first, checking the diagonal and the 2x2 condition
A=infsup(A,A);
for i=1:n
    if A.sup(i,i)<=0                                  % diagonal entry nonpositive
        pd=0;                                         % verified not positive definite
        x=zeros(n,1); x(i)=1;
        setround(gr); return
    end
    for j=i+1:n
        a=A(i,i)+A(j,j)-2*A(i,j);
        if a.sup<=0                                   % 2x2 condition not satisfied
            pd=0;                                     % verified not positive definite
            x=zeros(n,1); x(i)=1; x(j)=-1;
            setround(gr); return
        end  
        a=A(i,i)+A(j,j)+2*A(i,j);
        if a.sup<=0                                   % 2x2 condition not satisfied
            pd=0;                                     % verified not positive definite
            x=zeros(n,1); x(i)=1; x(j)=1;
            setround(gr); return
        end 
    end                                       
end
% second, finding a random x satisfying x'*A*x<=0
for i=1:max(n^2,1e03)
    x=2*rand(n,1)-1;
    while isequal(x,zeros(n,1))
        x=2*rand(n,1)-1;
    end
    a=x'*A*x;
    if a.sup<=0                                       % x'*A*x verified nonpositive
        pd=0;                                         % verified not positive definite
        setround(gr); return
    end  
end
x=repmat(NaN,n,1);                                    % no verified result
setround(gr);

   
