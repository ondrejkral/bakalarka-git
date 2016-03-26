function [x,y]=verlinineqnn(A,b)
%    VERLININEQNN       Verified nonnegative solution of a system of linear inequalities.  
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a rectangular real matrix A and a matching real vector b,
%    [x,y]=verlinineqnn(A,b)
%    either computes a verified solution of the system of linear inequalities
%        A*x <= b,       (1)
%          x >= 0,       (2)
%    or verifies nonexistence of a solution, or yields no verified result.
%
%    Possible outputs:
%
%    ~isnan(x(1)) :              x is a real vector verified to satisfy (1), (2), and 
%                                y is a vector of NaN's,
%    ~isnan(y(1)) :              y is a real vector verified to satisfy A'*y>=0, y>=0, b'*y<=-1 
%                                (which by Farkas lemma implies nonexistence of a solution to (1), (2)), and 
%                                x is a vector of NaN's,
%    isnan(x(1)) & isnan(y(1)) : no verified result.
%
%    See also LINPROG.

%    Copyright 2008 Jiri Rohn
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
%    2008-01-05   first version 
%    2008-01-16   version for posting
%
gr=getround;
setround(0);
b=b(:); [m,n]=size(A);  
x=repmat(NaN,n,1); 
y=repmat(NaN,m,1);
if nargin~=2||m~=length(b)||~isreal(A)||~isreal(b)||isintval(A)||isintval(b)
    setround(gr); return
end
if ~issparse(A)
    A=sparse(A);
end
xx=verlinineqnninner(A,b); 
if ~isnan(xx(1))
    x=xx; % verified solution of A*x<=b, x>=0
    setround(gr); return
end
Ao=[-A'; -speye(m,m); b'];
bo=[zeros(1,n+m) -1]';
yy=verlinineqnninner(Ao,bo);
if ~isnan(yy(1))
    y=yy; % verified solution of A'*y>=0, y>=0, b'*y<=-1 (implies nonexistence of solution to A*x<=b, x>=0)
    setround(gr); return
end
setround(gr);
%
% Subfunctions
%
function x=verlinineqnninner(A,b)
% inner subroutine of verlinineqnn
% finds a verified solution to A*x<=b, x>=0, or yields a vector of NaN's
% additive and multiplicative perturbation used
[m,n]=size(A);  
x=repmat(NaN,n,1); 
ep=max(1e-10,max([m n 100])*max([norm(A,inf) norm(b,inf)])*eps); 
e=ones(n,1);
Ao=[A; -speye(n,n)];
bo=[b' zeros(n,1)']'; % Ao*x<=bo is equivalent to A*x<=b, x>=0
% additive perturbation
bo=bo-ep*ones(m+n,1);
xx=lpprocedure(e,Ao,bo); % solves min e'*x subject to Ao*x<=bo
xi=infsup(xx,xx); % interval quantity
left=A*xi;
if all(left.sup<=b)&&all(xx>=0)
    x=xx; % real quantity; verified solution
    return
end
% multiplicative perturbation
for i=1:m+n
    if bo(i)~=0
        bo(i)=bo(i)-ep*abs(bo(i));
    else
        bo(i)=bo(i)-ep;
    end
end
xx=lpprocedure(e,Ao,bo); % solves min e'*x subject to Ao*x<=bo
xi=infsup(xx,xx); % interval quantity
left=A*xi;
if all(left.sup<=b)&&all(xx>=0) 
    x=xx; % real quantity; verified solution
    return
end
%
function x=lpprocedure(c,A,b)
% solves linear programming problem min c'*x subject to A*x<=b
% x should be always assigned (unverified optimal solution, or something else; the result is checked afterwards)
% placed separately so that a different linear programming procedure might also be used
wg=warning;                
warning off
options=optimset('linprog');                                  
options=optimset(options,'Display','off'); % suppresses displaying linprog's error messages
x=linprog(c,A,b,[],[],[],[],[],options);
warning(wg);


