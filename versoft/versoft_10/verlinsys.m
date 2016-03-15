function [sol,x,B,E]=verlinsys(A,b)
%    VERLINSYS       Verified description of all solutions of a system of linear equations.   
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For an m-by-n real matrix A and a matching real vector b,
%        [sol,x,B,E]=verlinsys(A,b)
%    yields a verified description of ALL solutions of the system of linear
%    equations A*x=b. The number of solutions is given in the variable sol:
%    
%    sol = Inf     A*x=b is verified to have infinitely many solutions; 
%                  x is an n-by-1 interval vector and B is an n-by-n
%                  interval matrix that are verified to contain a real
%                  vector xo and a real matrix Bo such that the set X of
%                  all solutions of A*x=b is described by 
%                      X = { xo+Bo*y | y in R^n }; 
%                  moreover, x has at least n-r zero entries, where r is
%                  the rank of A ("basic solution"),
%    sol =  1      A*x=b is verified to have exactly one solution xo which
%                  is verified to be contained in the interval vector x; 
%                  B is a matrix of interval zeros,
%    sol =  0      A*x=b is verified to possess no solution; x and B
%                  consist of NaN's, 
%    sol = -1      no verified result; x and B consist of NaN's.
%
%    The structure E explains reasons for NaN output. 
%
%    If sol==0, you may try VERLSQ for verified least squares solution(s).
%
%    See also VERPINV, VERBASIS, VERFULLCOLRANK, VERLSQ, VERIFYLSS.

%    Copyright 2008 Jiri Rohn
%
%    Based on the description of the set X of all solutions (if nonempty) as
%        X = { xo+(eye(n,n)-pinv(A)*A)*y | y in R^n },
%    where xo in X, with a verified pseudoinverse being computed by VERPINV.  
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
%    2008-04-14   first version 
%    2008-04-20   version for posting
%
gr=getround;
setround(0);
b=b(:); 
[m,n]=size(A); 
% defaults
sol=-1;
x=repmat(infsup(NaN,NaN),n,1);
B=repmat(infsup(NaN,NaN),n,n);
E.error='verlinsys: none';
E.where='NaN';
E.value='NaN';
% data check
if nargin~=2||nargout>4||m~=length(b)||~isreal(A)||~isreal(b)||isintval(A)||isintval(b)
    E.error='verlinsys: wrong data';
    setround(gr); return
end
% square case first
if m==n
    xx=verifylss(A,b);
    if ~isnan(xx.inf(1)) % computed
        sol=1;
        x=xx;
        B=repmat(infsup(0,0),n,n); % B set to zero interval matrix
        setround(gr); return
    end
    % not computed
end
% quantities I (for sol==0 and sol==1)
X=verpinv(A); % pseudoinverse
if isnan(X.inf(1,1))
    E.error='verlinsys: verified pseudoinverse of A could not be computed';
    setround(gr); return % X needed in all three conditions
end
res=(eye(m,m)-A*X)*b; % should be ~=0 for nonexistence of solution
BB=eye(n,n)-X*A; % should be nonzero for infinitely many solutions
fcrA=verfullcolrank(A); % should be fcrA==1, fcrAb==0 for existence of solution
fcrAb=verfullcolrank([A b]);
% no solution
if ~isnan(X.inf(1,1))&&(any(res.sup<0)||any(res.inf>0)) % verified no solution
    sol=0;
    setround(gr); return
end
% unique solution
if ~isnan(X.inf(1,1))&&fcrA==1&&fcrAb==0 % verified unique solution
    sol=1;
    x=X*b;
    B=repmat(infsup(0,0),n,n); % B set to zero interval matrix
    setround(gr); return
end
% quantities II (for sol==Inf)
[BS,K]=verbasis(A); % basis 
if isnan(BS(1,1)) 
    E.error='verlinsys: verified pseudoinverse of the basis could not be computed';
    setround(gr); return 
end
fcrBSb=verfullcolrank([BS b]); % should be fcrBSb==0 for existence of solution
XBS=verpinv(BS); % used for computation of x
% infinitely many solutions
if ~isnan(X.inf(1,1))&&~isnan(BS(1,1))&&~isnan(XBS.inf(1,1))&&fcrBSb==0&&(any(any(BB.inf>0))||any(any(BB.sup<0))) 
    sol=Inf;
    xK=XBS*b; % xK verified to enclose a solution of BS*X=b
    x=infsup(zeros(n,1),zeros(n,1)); x(K)=xK; % x verified to enclose a solution of A*X=b
    B=BB;
    setround(gr); return
end
E.error='verlinsys: quantities needed could not be computed; no verified result';
setround(gr);

