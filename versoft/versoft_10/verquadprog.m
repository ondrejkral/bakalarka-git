function [flag,x,y,h]=verquadprog(A,b,c,D)
%    VERQUADPROG     Verified convex quadratic programming.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For real data A, b, c, D, with D positive definite (A, D full or sparse),
%    [flag,x,y,h]=verquadprog(A,b,c,D)
%    either computes a verified optimal solution x, verified Lagrangian
%    multipliers y and a verified optimal value h of the quadratic
%    programming problem 
%        min 0.5*x'*D*x+c'*x   subject to   A*x>=b, x>=0, 
%    or verifies (in)feasibility, or yields no verified result. The
%    respective outcome is always described verbally in the variable flag. 
%
%    Possible outputs:
%
%    flag='verified optimum   ' : x is verified to enclose an optimal solution, 
%                                 y is verified to enclose a vector of Lagrangian multipliers, 
%                                 h is verified to enclose the optimal value,
%    flag='verified feasible  ' : x is verified to enclose a feasible solution
%                                   (optimality nor unboundedness could be verified),
%                                 y, h are NaN's,
%    flag='verified infeasible' : y is a Farkas vector verified to satisfy y>=0, A'*y>=0, b'*y<0
%                                   (whose existence proves infeasibility),
%                                 x, h are NaN's,
%    flag='verified not PD    ' : D verified not to be positive definite,    
%                                   x is a vector verified to satisfy x'*D*x<=0, 
%                                   y, h are NaN's,
%    flag='sizes do not match ' : x, y, h are NaN's
%                                   (sizes of A, b, c, D do not match, or D is nonsymmetric),
%    flag='no verified result ' : x, y, h are NaN's
%                                   (no verified result could be found).
%
%    Complexity: The algorithm solves one quadratic programming problem and
%    at most three linear programming problems (independently of the size
%    of the original problem), and employs two verification procedures
%    which both run approximately in O(m^3) time, where m=size(A,1).  
%
%    See also QUADPROG, LINPROG, VERLINPROG, VERIFYLSS.

%    Copyright 2007 Jiri Rohn
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
warning('off');
% checking data
b=b(:); c=c(:); [m,n]=size(A);
flag='no verified result ';
x=repmat(infsup(NaN,NaN),n,1); y=repmat(infsup(NaN,NaN),m,1); h=infsup(NaN,NaN);
if ~(length(b)==m&&length(c)==n&&isequal(size(D),[n,n])&&isequal(D,D')) 
    flag='sizes do not match ';
    setround(gr); return                                      
end
if isintval(A)||isintval(b)||isintval(c)||isintval(D)           % data not thin
    setround(gr); return
end
% verifying positive definiteness
[pd,xpd]=verpd(D);
if pd==0                                                        % D verified not positive definite
    flag='verified not PD    '; x=xpd;
    setround(gr); return
end
if pd==-1||pd==-2
    setround(gr); return
end
Ao=[A -eye(m,m)];                                               
% verifying infeasibility
yi=verinfeas(Ao,b);
if ~isnan(yi.inf(1))                                            % verified Farkas vector found
    y=yi; flag='verified infeasible';
    setround(gr); return
end
% verifying feasibility
xf=veropt(Ao,b,ones(n+m,1));
if isnan(xf.inf(1))                                             % verified feasible solution not found
    flag='no verified result ';
    setround(gr); return
end
xf=xf(1:n);
% computing unverified optimum
options=optimset('quadprog');                                  
options=optimset(options,'Display','off');
[x,h,exitflag,output,lambda]=quadprog(D,c,-A,-b,[],[],zeros(n,1),Inf*ones(n,1),[],options);
    % solves min 0.5*x'*D*x + c'*x  subject to  A*x>=b, x>=0
if exitflag<=0
    x=xf; flag='verified feasible  ';                           % previous feasible solution outputed
    setround(gr); return
end
p=lambda.ineqlin;
% verifying optimality
M=[D -A'; A zeros(m,m)];                                        % data of the respective LCP
q=[c' -b']';
z=[x' p']';                                                     % y, z computed in this way satisfy (in infinite precision)
y=M*z+q;                                                        % y=M*z+q, y>=0, z>=0, y'*z=0
x=y-z;
I=eye(m+n,m+n);
I=infsup(I,I);                                                  % x satisfies (in infinite precision)
M=infsup(M,M);                                                  % (I+M)*x+(I-M)*abs(x)=2*q
A=I+M; B=I-M; b=infsup(2,2)*infsup(q,q);
xx=absvalverifn(A,B,b,x);                                       
if isnan(xx.inf(1))
    x=xf; flag='verified feasible  ';                           % previous feasible solution outputed
    setround(gr); return                                        % no verified solution of LCP found
end
% extracting x^-=max(-x,0)                                      % xx is verified solution
for i=1:m+n                                        
    if 0<xx.inf(i)
        xx(i)=infsup(0,0);
    end
    if xx.sup(i)<0
        xx(i)=-xx(i);
    end
    if (xx.inf(i)<=0)&&(xx.sup(i)>=0)
        xx(i)=infsup(0,-xx.inf(i));
    end
end
% verified quantities
x=xx(1:n); p=xx(n+1:n+m); y=p;                                  % verified optimal solution and Lagrangian multipliers
D=infsup(D,D); c=infsup(c,c);              
h=0.5*(x'*D*x)+c'*x;                                            % verified optimal value
flag='verified optimum   ';
setround(gr);
%
% Subfunctions
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
if isspd(A)==1                                                  % routine isspd by S. M. Rump
    pd=1;                                                       % verified positive definite
    setround(gr); return
end
% verifying not-positive-definiteness
% first, checking the diagonal and the 2x2 condition
A=infsup(A,A);
for i=1:n
    if A.sup(i,i)<=0                                            % diagonal entry nonpositive
        pd=0;                                                   % verified not positive definite
        x=zeros(n,1); x(i)=1;
        setround(gr); return
    end
    for j=i+1:n
        a=A(i,i)+A(j,j)-2*A(i,j);
        if a.sup<=0                                             % 2x2 condition not satisfied
            pd=0;                                               % verified not positive definite
            x=zeros(n,1); x(i)=1; x(j)=-1;
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
    if a.sup<=0                                                 % x'*A*x verified nonpositive
        pd=0;                                                   % verified not positive definite
        setround(gr); return
    end  
end
x=repmat(NaN,n,1);                                              % no verified result
setround(gr);
%
function [x,B,N]=veropt(A,b,c)                                  % unchanged against verlinprog
% B is the "basis index set" of an optimal solution of the LP problem
% min c'*x  subject to  A*x=b, x>=0,
% x is a verified basic feasible solution with this basis
% N is the set of nonbasic indices
[m,n]=size(A); 
x=repmat(infsup(NaN,NaN),n,1);
B=repmat(NaN,m,1); N=repmat(NaN,n,1);
options=optimset('linprog');                                  
options=optimset(options,'Display','off');                     % suppresses displaying linprog's error messages
[xopt,fval,exitflag]=linprog(c,-eye(n),zeros(n,1),A,b,[],[],[],options);
if exitflag<=0 
    return
end
[xx,J]=sort(xopt); B=J(n-m+1:n); N=J(1:n-m);                    % B is set of "basic" indices, N of "nonbasic" ones
% xB=verifylss(A(:,B),b);                                       % old version using full A; replaced by next 5 lines
AB=A(:,B);
if issparse(AB)
    AB=full(AB);                                                % only the square submatrix taken full
end
xB=verifylss(AB,b);
if isnan(xB.inf(1))||~all(xB.inf>=0),                           % verified "optimal" solution not found
    return
end
x=infsup(zeros(n,1),zeros(n,1)); x(B)=xB;                       % verified "optimal" solution found
%
function y=verinfeas(A,b)                                       % unchanged against verlinprog
% y is a Farkas vector (i.e., satisfying A'*y>=0, b'*y<0)
% its existence implies infeasibility of A*x=b
[m,n]=size(A);
y=repmat(infsup(NaN,NaN),m,1);
ep=max(1e-08,max([m n 100])*max([norm(A,inf) norm(b,inf)])*eps);
Afv=[A' -A'    -eye(n) zeros(n,1);                              % Afv is (n+1)x(2*m+n+1)
     b' -b' zeros(1,n)         1]; 
bfv=[zeros(n,1)' -1]';                                          % bfv is (n+1)x1
bfv=bfv+ep*[ones(1,n) -1]';                                     % perturbation to compensate roundoff errors (so that A'*y>=0)
yf=veropt(Afv,bfv,ones(2*m+n+1,1));                             % system: A'*y>=0, b'*y<=-1, y written as y=y1-y2
if ~isnan(yf.inf(1))
    yf=mid(yf);
    y1=yf(1:m);
    y2=yf(m+1:2*m);
    yf=y1-y2;                                                   % would-be Farkas vector
    A=infsup(A,A); b=infsup(b,b); yf=infsup(yf,yf);             % (i.e., should satisfy A'*y>=0, b'*y<0)
    alpha=A'*yf; beta=b'*yf;
    if (all(alpha.inf>=0))&&(beta.sup<0)                        % infeasibility verified
        y=yf;                                                   % Farkas vector outputed
        return
    else
        return
    end
else
    return
end
%
function xx=absvalverifn(A,B,b,x)                     % absvalverifn of 2007-03-26
%    ABSVALVERIFN     Verification of a solution of the equation A*x+B*abs(x)=b.
%
%    For an approximate real solution x of the equation
%        A*x+B*abs(x)=b,     (1)
%    xx=absvalverifn(A,B,b,x)
%    either produces a tight interval vector xx verified to contain a solution 
%    of (1), or fails (yields an interval vector xx of NaN's).

%    Copyright 2007 Jiri Rohn
%
n=size(A,1);
xx=repmat(infsup(NaN,NaN),n,1);
if ~isintval(A), A=infsup(A,A); end                   % converting data from real to intval 
if ~isintval(B), B=infsup(B,B); end
if ~isintval(b), b=infsup(b,b); end
if  isintval(x), return, end
if issparse(A)                                        % A, B made full because of verifylss
    A=full(A);                                               
end
if issparse(B)
    B=full(B);                                               
end
I=eye(n,n);
z=ones(n,1);
for j=1:n
    if x(j)<0, z(j)=-1; end                           % sign vector of x
end
infinf=infsup(repmat(-Inf,n,1),repmat(Inf,n,1));
x1=verifylss(A+B*diag(z),b);                          % enclosure via verifylss
if ~(~isnan(x1.inf(1))&&all(z.*inf(x1)>=0)&&all(z.*sup(x1)>=0))
    x1=infinf;                                        % failure to produce verified output
end
M=inv(I-abs(inv(A)*B));
if (~isnan(M.inf(1,1)))&&(all(all(M.inf>=0)))         % M verified nonnegative
    rad=M*abs(inv(A)*(A*x+B*abs(x)-b));
    x2=infsup(x-rad.sup,x+rad.sup);                   % enclosure via |x^*-x|<=M*|inv(A)*residual|
else
    x2=infinf;                                        % failure to produce verified output
end
x3=intersect(x1,x2);                                  % intersection of the two enclosures
if ~isequal(x3,infinf)                                % at least one enclosure computed
    xx=x3;                                            % verified enclosure of the solution of Ax+B|x|=b
end
