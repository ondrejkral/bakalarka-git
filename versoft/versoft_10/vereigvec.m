function [evc,lambda,As]=vereigvec(A,x)
%    VEREIGVEC      Verified real eigenvector of an interval matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a square interval matrix A and a REAL vector x,
%    [evc,lambda,As]=vereigvec(A,x)
%    verifies x to be an eigenvector of some matrix in A,
%    or not to be an eigenvector of any matrix in A,
%    or yields no verified result (unfortunately, complex eigenvectors
%    cannot be handled yet): 
%
%    evc= 1           x is verified to be an eigenvector of some matrix in A,
%                     lambda is an interval number such that for each
%                       lambda0 in lambda, A is verified to contain a matrix
%                       having (lamda0,x) as an eigenpair,
%                     As is a very tight interval matrix verified to contain 
%                       a matrix having (mid(lambda),x) as an eigenpair,
%    evc= 0           x is verified not to be an eigenvector of any matrix in A,
%                     lambda and As consist of NaN's, 
%    evc=-1           no verified result (data may be wrong).
%
%    See also VEREIGVAL.

%    Copyright 2007 Jiri Rohn
%
%    Based on the section "Real eigenvectors" in
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
% checking data
x=x(:);
[m,n]=size(A); p=length(x);
evc=-1; lambda=infsup(NaN,NaN); As=repmat(infsup(NaN,NaN),n,n);
if m~=n||n~=p||~isreal(A)||~isreal(x)||isequal(x,zeros(p,1)) % error('wrong data')
   setround(gr); return
end
if ~isintval(A), 
    A=infsup(A,A);                                    % allows for real input 
end
% checking the basic inequality
[ac,Delta]=vermidrad(A);
z=sgn(x);
x1=infsup(x,x);                                       % x double, x1 intval
Tz=diag(z); 
left =Tz*(ac-Tz*Delta*Tz)*x1*x1'*Tz;                  % left-hand  side of the inequality
right=Tz*x1*x1'*(ac+Tz*Delta*Tz)'*Tz;                 % right-hand side of the inequality
% inequality verified not to be satisfied
if any(any(right.sup<left.inf))                       % verified not to be an eigenvector
    evc=0; 
    setround(gr); return
end
% inequality verified to be satisfied
if all(all(left.sup<=right.inf))                      % verified to be an eigenvector; Rohn, SIMAX 1993, Thm. 4.1
    B=find(x~=0);
    denleft= (Tz*ac*Tz-Delta)*abs(x);                 
    denright=(Tz*ac*Tz+Delta)*abs(x);                
    num=abs(x);                                       
    denleft=denleft(B);
    denright=denright(B);
    num=num(B);
    left= denleft./num;                               % left  ratio
    right=denright./num;                              % right ratio
    lambdal=max(left.sup);                            % verified lower bound of lambda        
    lambdau=min(right.inf);                           % verified upper bound of lambda 
    if lambdal>lambdau                                % bounds contradict: no verified solution
        evc=-1; 
        lambda=infsup(NaN,NaN); As=repmat(infsup(NaN,NaN),n,n); 
        setround(gr); return
    end
    lambda=infsup(lambdal,lambdau);                   % lambda
    lambdam=mid(lambda);                              % midpoint of lambda
    % finding a matrix with eigenpair (lambdam,x)
    A1=A;
    for i=1:n
        A1(i,i)=A1(i,i)-lambdam;                      % A1=A-lambdam*I
    end
    AAs=versingnull(A1,x);                            % enclosure of a singular matrix in A1 having x as a null vector
    if isnan(AAs.inf(1,1))                            % no enclosure outputed
        evc=-1;                                       % no verified result
        lambda=infsup(NaN,NaN); As=repmat(infsup(NaN,NaN),n,n);          
        setround(gr); return
    end
    for i=1:n
        AAs(i,i)=AAs(i,i)+lambdam;                    % AAs=AAs+lambdam*I (back to A)
    end
    if in(AAs,A)                                      % AAs part of A
        evc=1; As=AAs;                                % (lambdam,x) is an eigenpair of a matrix in As; Rohn, SIMAX 1993, proof of Thm. 4.1
        setround(gr); return                          
    else
        evc=-1;                                       % AAs not a part of A: no verified result
        lambda=infsup(NaN,NaN); As=repmat(infsup(NaN,NaN),n,n);          
        setround(gr); return
    end
end
evc=-1;                                               % no verified result
lambda=infsup(NaN,NaN); As=repmat(infsup(NaN,NaN),n,n); 
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
%
function As=versingnull(A,x)
%    VERSINGNULL    Verified singular matrix in A having x as a null vector.
% 
%    ~isnan(As.inf(1,1)): As is a tight interval matrix verified to be a part
%                         of A and to contain a singular matrix having x
%                         as a null vector
%     isnan(As.inf(1,1)): no result
%
[m,n]=size(A);
As=repmat(infsup(NaN,NaN),n,n);
if m~=n
    return
end
if nargin==1||isnan(x(1))
    return
end
z=sgn(x);
xi=infsup(x,x);
[Ac,Delta]=vermidrad(A);
oeprl=abs(Ac*xi);                                     % Oettli-Prager inequality, left  side
oeprr=Delta*abs(xi);                                  % Oettli-Prager inequality, right side
if all(oeprl.sup<=oeprr.inf)                          % Oettli-Prager inequality satisfied, singularity of A verified
    y=(Ac*xi)./oeprr;
    y(find(isnan(y)))=infsup(1,1);                    % case od both numerator and denominator being zero
    As=Ac-(diag(y)*Delta)*diag(z);                    % construction of singular As ...
    As=intersect(As,A);                               % ... in A
    if ~any(any(isnan(As)))                           % intersection nowhere empty
        return                                        % with output As
    else
        As=repmat(infsup(NaN,NaN),n,n); 
        return                                        % with As of NaN's, but still verified singular (this fact not used here)
    end
else
    return
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

