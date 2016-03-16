function [x,As]=verintlinineqs(A,b)
%    VERINTLININEQS     Verified strong solution of interval linear inequalities. 
%
%    For a rectangular interval matrix A and a matching interval vector b, 
%    [x,As]=verintlinineqs(A,b)
%    either computes a strong solution x to A*X<=b (i.e., a real vector x
%    verified to satisfy Ao*x<=bo for each Ao in A and bo in b), or verifies
%    nonexistence of such a solution, or yields no verified result.
%
%    Possible outputs:
%
%    ~isnan(x(1)) :                      x is a verified strong solution of A*X<=b,
%                                        and As is an interval matrix of NaN's,
%    ~isnan(As.inf(1,1)) :               As is a very tight ("almost thin") interval 
%                                        matrix verified to contain a real matrix Ao 
%                                        such that the system Ao*x<=b.inf has no solution 
%                                        (which proves that no strong solution exists), 
%                                        and x is a vector of NaN's,
%     isnan(x(1)) & isnan(As.inf(1,1)) : no verified output.
%
%    COMMENT. A theoretical result [1] asserts that if each system Ao*x<=bo, 
%    where Ao in A and bo in b, has a solution (depending generally on Ao and bo), 
%    then there exists a vector x satisfying Ao*x<=bo for EACH Ao in A and bo in b. 
%    Such a vector x is called a strong solution of the system A*X<=b. 
%
%    [1] J. Rohn and J. Kreslova, Linear Interval Inequalities, LAMA 38 (1994), 79-82.
%
%    See also VERLININEQNN.

%    Copyright 2008 Jiri Rohn
%
%    Based on Section 2.13 in 
%    M. Fiedler, J. Nedoma, J. Ramik, J. Rohn and K. Zimmermann, Linear Optimization 
%    Problems with Inexact Data, Springer-Verlag, New York 2006.
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
%    2008-01-20   version for posting
%
gr=getround;
setround(0);
b=b(:); [m,n]=size(A);  
x=repmat(NaN,n,1); As=repmat(infsup(NaN,NaN),m,n);
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
Al=inf(A); Au=sup(A); bl=inf(b); % the bounds
Ao=[Au -Al]; % matrix of the system; see Fiedler et al., (2.89)
[xx,y]=verlinineqnn(Ao,bl); % finds verified nonnegative solution of Ao*x<=bl
if ~isnan(xx(1)) % solution found
    xxi=infsup(xx,xx);
    xxi=xxi(1:n)-xxi(n+1:2*n); % interval vector of the original size
    X=[xx(1:n)-xx(n+1:2*n) xxi.inf xxi.mid xxi.sup]; % noninterval vectors; candidates for strong solution
    [Ac,Delta]=vermidrad(A);
    [bc,delta]=vermidrad(b);
    for x1=X
        left =Ac*x1-bc;
        right=-Delta*abs(x1)-delta;
        if all(left.sup<=right.inf) % Fiedler et al., (2.94); strong solution found
            x=x1; % verified strong solution
            setround(gr); return 
        end
    end
    setround(gr); return % no result 
end
if ~isnan(y(1)) % Ao*x<=bl verified not to have a nonnegative solution
    As=vernull(A',y); % Fiedler et al., proof of Thm. 2.23
    if ~isnan(As.inf(1,1)) 
        As=full(As'); % Ao*x<=bl unsolvable for some Ao in As which is a part of A
        setround(gr); return
    end
    setround(gr); return
end
setround(gr); % no result
%
% Subfunctions
%
function As=vernull(A,x)
%    VERNULL    Verified matrix in A having x as a null vector.
%
%    ~isnan(As.inf(1,1)): As is a tight interval matrix verified to be a part of A 
%                         and to contain a thin matrix Ao having x as a null vector
%                         (i.e., Ao*x=0),
%     isnan(As.inf(1,1)): no result.
%
[m,n]=size(A); p=length(x);
As=repmat(infsup(NaN,NaN),m,n);
if n~=p||nargin~=2||any(isnan(x))||~isintval(A)||isintval(x)
    return
end
z=sgn(x);
xi=infsup(x,x);
[Ac,Delta]=vermidrad(A);
oeprl=abs(Ac*xi);                                     % Oettli-Prager inequality, left  side
oeprr=Delta*abs(xi);                                  % Oettli-Prager inequality, right side
if all(oeprl.sup<=oeprr.inf)                          % Oettli-Prager inequality satisfied, x verified null vector of A
    y=(Ac*xi)./oeprr;
    y(find(isnan(y)))=infsup(1,1);                    % case od both numerator and denominator being zero
    As=Ac-(diag(y)*Delta)*diag(z);                    % construction of As ...
    As=intersect(As,A);                               % ... in A
    if ~any(any(isnan(As)))                           % intersection nowhere empty
        return                                        % with output As
    else
        As=repmat(infsup(NaN,NaN),m,n); 
        return                                        % with As of NaN's, but x still verified null vector of A
    end
else
    return
end
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
%
function z=sgn(x)
% signum of x for real or intval x
% x column or row, z column
n=length(x);
z=zeros(n,1);
if ~isintval(x) % real vector
    for i=1:n
        if x(i)>=0
            z(i)=1;
        else
            z(i)=-1;
        end
    end
else % intval vector
    for i=1:n
        if x.inf(i)>0
            z(i)=1;
        end
        if x.sup(i)<0
            z(i)=-1;
        end
        % otherwise z(i)=0
    end
end




