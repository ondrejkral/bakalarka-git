function [pd,x,E]=verpd(A)
%    VERPD    Verified positive definiteness of a real matrix.
%
%    For a symmetric real (not interval) matrix A,
%        [pd,x,E]=verpd(A)
%    verifies positive definiteness or not-positive-definiteness of A,
%    or yields no verified result:
%
%    pd= 1 :     A verified positive definite,
%    pd= 0 :     A verified not to be positive definite, and x is a nonzero real
%                vector verified to satisfy x'*A*x<=0 (a "certificate"), 
%    pd=-1 :     no verified result.
%
%    If pd~=0, then x is a vector of NaN's.
%
%    The structured array E explains reasons for NaN output. It has three fields: 
%    E.error, E.where, E.value. 
%
%    See also ISSPD, VERPOSDEF.

%    Copyright 2007-2008 Jiri Rohn.
%    The routine ISSPD copyrighted by Siegfried M. Rump.
%
%    Positive definiteness is first checked by ISSPD. If the test fails,
%    not-positive-definiteness is checked by repeated use of the 2-by-2
%    condition. If this fails, VEREIG is brought in to judge
%    (not-)positive definiteness from verified eigenvalues. If even this
%    fails, finally an attempt is made at finding through random search
%    a real vector x verified to satisfy x'*A*x<=0.
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
%    History
%
%    2007, Spring   created as a part of VERQUADPROG
%    2008-05-20     made independent
%    2008-05-21     eigenvalue part added
%    2008-11-16     output variable E added, vereig called instead of ol
%    2008-11-22     final version (html copy)
%
gr=getround;
setround(0);
% checking data
n=size(A,1);
pd=-1; 
x=repmat(NaN,n,1);
E.error='verpd: none';
E.where='NaN';
E.value='NaN';
if ~isreal(A)||isintval(A)
    E.error='verpd: matrix not real or of type intval';
    setround(gr); return
end
if ~isequal(A,A')
    E.error='verpd: matrix not symmetric';
    setround(gr); return
end
% verifying positive definiteness via isspd
if isspd(A)==1                                        % routine isspd by S. M. Rump
    pd=1;                                             % verified positive definite
    setround(gr); return
end
% verifying not-positive-definiteness via the 2x2 condition (see e.g. Golub and van Loan)
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
% verifying (not-)positive-definiteness via vereig
[L,X]=vereig(A);
lam=diag(L);
if ~isnan(lam.inf(1))                                 % eigenvalues computed
    if lam.inf(1)>0
        pd=1;                                         % verified positive definite (all eigenvalues verified positive)
        setround(gr); return
    end
    if lam.sup(1)<=0                                  
        a=X(:,1)'*A*X(:,1);
        if a.sup<=0                                   
            pd=0;                                     % verified not positive definite (minimal eigenvalue verified nonpositive)
            x=X(:,1);
            setround(gr); return
        end 
    end
end
% finally, trying to find a random x satisfying x'*A*x<=0
for i=1:max(n,1e03)
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
% all tests failed
pd=-1;
x=repmat(NaN,n,1);                                    % no verified result
E.error='verpd: failure';
setround(gr);  
% end of verpd
%
