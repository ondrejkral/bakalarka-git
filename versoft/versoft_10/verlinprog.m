function [flag,x,y,h]=verlinprog(A,b,c)
%    VERLINPROG       Verified linear programming.  
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a real matrix A (full or sparse) and matching real vectors b, c,
%    [flag,x,y,h]=verlinprog(A,b,c)
%    either computes verified optimal solution x, verified dual optimal solution y
%    and verified optimal value h of the linear programming problem
%        min c'*x   subject to   A*x=b, x>=0,
%    or verifies (in)feasibility, or verifies unboundedness, or yields no
%    verified result. The respective outcome is always described verbally
%    in the variable flag.
%
%    Possible outputs:
%
%    flag='verified optimum   ' : x is verified to enclose a primal optimal solution, 
%                                 y is verified to enclose a dual optimal solution, 
%                                 h is verified to enclose the optimal value,
%    flag='verified unbounded ' : x is verified to enclose a primal feasible solution xo, and
%                                 y is verified to enclose a vector yo such that the objective 
%                                   tends to -Inf along the feasible half-line {xo+t*yo | t>=0},
%                                 h is NaN,
%    flag='verified feasible  ' : x is verified to enclose a primal feasible solution
%                                   (optimality nor unboundedness could be verified),
%                                 y, h are NaN's,
%    flag='verified infeasible' : y is verified to enclose a Farkas vector yo satisfying A'*yo>=0, b'*yo<0
%                                   (whose existence proves primal infeasibility),
%                                 x, h are NaN's,
%    flag='no verified result ' : x, y, h are NaN's
%                                   (no verified result could be found).
%    flag='sizes do not match ' : x, y, h are NaN's
%                                   (sizes of A, b, c do not match).
%
%    Complexity: The algorithm solves at most four linear programming
%    problems (independently of the size of the original problem) and uses 
%    a verification procedure which runs approximately in O(m^3) time, where 
%    m=size(A,1).  
%
%    See also LINPROG, VERIFYLSS.

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
b=b(:); c=c(:); [m,n]=size(A); p=length(b); q=length(c);
flag='no verified result ';
x=repmat(infsup(NaN,NaN),n,1); y=repmat(infsup(NaN,NaN),m,1); h=infsup(NaN,NaN); 
if ~(m==p&&n==q)||(m>n), 
    flag='sizes do not match ';
    setround(gr); return
end
if isintval(A)||isintval(b)||isintval(c),                       % error('data not real')
    setround(gr); return
end
if issparse(b)
    b=full(b);
end
if issparse(c)
    c=full(c);
end
% verifying infeasibility
yi=verinfeas(A,b);
if ~isnan(yi.inf(1))                                            % verified Farkas vector found
    y=yi; flag='verified infeasible';
    setround(gr); return
end
% verifying feasibility
xf=veropt(A,b,ones(n,1));
if isnan(xf.inf(1))                                             % verified feasible solution not found
    flag='no verified result ';
    setround(gr); return
end
% verifying unboundedness
yu=verunbound(A,c);
if ~isnan(yu.inf(1))                                            % verified descent direction found
    x=xf; y=yu; flag='verified unbounded ';
    setround(gr); return
end
% verifying optimality
[xo,B,N]=veropt(A,b,c);
if isnan(xo.inf(1))                                             % verified feasible primal solution with basis B not found
    x=xf; flag='verified feasible  ';                           % previous feasible solution outputed
    setround(gr); return
end
% yB=verifylss(A(:,B)',c(B));                                   % old version using full A; replaced by next 5 lines
AB=A(:,B);
if issparse(AB)
    AB=full(AB);                                                % only the square submatrix taken full
end
yB=verifylss(AB',c(B));
if isnan(yB.inf(1))                                             % verified feasible dual solution not found
    x=xo; flag='verified feasible  ';                           % candidate for optimum outputed as feasible solution
    setround(gr); return
end                                
c=infsup(c,c); A=infsup(A,A);
crit=c'-yB'*A;                                                  % criterial row (dual feasibility)
crit=crit(N);                                                   % nonbasic part of it
if ~all(crit.inf>=0)                                            % verified feasible dual solution not found
    x=xo; flag='verified feasible  ';                           % candidate for optimum outputed as feasible solution
    setround(gr); return
end                               
% verified quantities                                           % verified primal and dual feasible solutions found
x=xo;                                                           % x is a verified primal optimal solution
y=yB;                                                           % y is a verified dual optimal solution
h1=c'*x; h2=b'*y;                                               
h=intersect(h1,h2);                                             % h is a verified optimal value (duality theorem) 
if isnan(h.inf)
    h=h1;
end
flag='verified optimum   ';
setround(gr); 
%
% Subfunctions
%
function [x,B,N]=veropt(A,b,c)
% B is the "basis index set" of an optimal solution of the LP problem
% min c'*x  subject to  A*x=b, x>=0,
% x is a verified basic feasible solution with this basis
% N is the set of nonbasic indices
[m,n]=size(A); 
x=repmat(infsup(NaN,NaN),n,1);
B=repmat(NaN,m,1); N=repmat(NaN,n,1);
options=optimset('linprog');                                  
options=optimset(options,'Display','off');                      % suppresses displaying linprog's error messages
[xopt,fval,exitflag]=linprog(c,-eye(n),zeros(n,1),A,b,[],[],[],options);
if exitflag<=0 
    return
end
[xx,J]=sort(xopt); B=J(n-m+1:n); N=J(1:n-m);                    % B is set of "basic" indices, N of "nonbasic" ones
% xB=verifylss(A(:,B),b);                                       % old version using full A; replaced by next 5 lines
AB=A(:,B);
if issparse(AB)
    AB=full(AB);                                                % only the square submatrix taken full (because of verifylss)
end
xB=verifylss(AB,b);
if isnan(xB.inf(1))||~all(xB.inf>=0),                           % verified "optimal" solution not found
    return
end
x=infsup(zeros(n,1),zeros(n,1)); x(B)=xB;                       % verified "optimal" solution found
%
function y=verinfeas(A,b)
% y verified to enclose a Farkas vector yo (i.e., satisfying A'*yo>=0, b'*yo<0)
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
function y=verunbound(A,c)
% y verified to enclose a vector yo satisfying A*yo=0, yo>=0, c'*yo<=-1
% under feasibility its existence implies unboundedness
[m,n]=size(A);
y=repmat(infsup(NaN,NaN),n,1);
Aunb=[A zeros(m,1);                                             % Aunb is (m+1)x(n+1)
      c'        1];
bunb=[zeros(1,m) -1]';                                          % bunb is (m+1)x1
yunb=veropt(Aunb,bunb,ones(n+1,1));                             % yunb is (n+1)x1
if ~isnan(yunb.inf(1))
    y=yunb(1:n);                                                % y satisfies A*y=0, y>=0, c'*y=-1 
    return
end

