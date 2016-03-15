function [ts,x]=vertolsol(A,b)
%    VERTOLSOL     Verified tolerance solution of a system of interval linear equations. 
%
%    For a rectangular interval matrix A and a matching interval vector b,
%    [ts,x]=vertolsol(A,b)
%    either computes a verified tolerance solution to A*X=b (i.e., a real vector x
%    verified to satisfy Ao*x in b for each Ao in A), or verifies
%    nonexistence of such a solution, or yields no verified result.
%
%    Possible outputs:
%
%    ts= 1 :    x is a verified tolerance solution of A*X=b,
%    ts= 0 :    A*X=b is verified not to possess a tolerance solution,
%               and x is a vector of NaN's,
%    ts=-1 :    no verified output.
%
%    See also VERLININEQNN.

%    Copyright 2008 Jiri Rohn
%
%    Based on the section "Tolerance solutions" in
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
%    History
%
%    2007-02-22   first version 
%    2008-01-19   version for posting
%
gr=getround;
setround(0);
b=b(:); [m,n]=size(A);  
ts=-1; x=repmat(NaN,n,1); 
if nargin~=2||m~=length(b)||~isreal(A)||~isreal(b) 
    setround(gr); return
end
if ~isintval(A) % allows for real input
    A=infsup(A,A); 
end
if ~isintval(b) 
    b=infsup(b,b); 
end
if ~issparse(A)
    A=sparse(A); % makes A sparse, because of verlinineqnn
end
Al=inf(A); Au=sup(A); % the bounds
bl=inf(b); bu=sup(b);
Ao=[ Au  -Al;    % data for the system Ao*x<=bo; see Handbook, p. 45, (5.2)
    -Al   Au];
bo=[ bu' -bl']';
[xx,y]=verlinineqnn(Ao,bo); % finds verified nonnegative solution of Ao*x<=bo
if ~isnan(xx(1)) % solution found
    xxi=infsup(xx,xx);
    xxi=xxi(1:n)-xxi(n+1:2*n); % interval vector of the original size
    X=[xx(1:n)-xx(n+1:2*n) xxi.inf xxi.mid xxi.sup]; % noninterval vectors; candidates for tolerance solution
    for x1=X
        if all(in(A*x1,b)) % tolerance solution found
            ts=1; x=x1; % verified tolerance solution
            setround(gr); return 
        end
    end
    setround(gr); return % no result 
end
if ~isnan(y(1)) % Ao*x<=bo verified not to have a nonnegative solution
    ts=0; % A*X=b verified not to have a tolerance solution
    setround(gr); return
end
setround(gr); % no result

