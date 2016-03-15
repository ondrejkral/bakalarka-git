function [x,B,E]=verlsq(A,b)
%    VERLSQ       Verified description of all linear squares solutions of a system of linear equations.   
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For an m-by-n real matrix A and a matching real vector b,
%        [x,B,E]=verlsq(A,b)
%    finds an n-by-1 interval vector x and an n-by-n interval matrix B that are 
%    verified to contain a real vector xo and a real matrix Bo such that
%    the set X of ALL least squares solutions of the equation A*x=b is
%    described by 
%        X = { xo+Bo*y | y in R^n }.
%    In particular, xo is the minimum 2-norm least squares solution. If A
%    is verified to have full column rank (so that xo is the unique
%    solution), then B is the interval matrix of exact zeros. If no
%    verified output is given, then x, B consist of NaN's. 
%
%    The structure E explains reasons for NaN output. It has three fields:
%    E.error, E.where, E.value.  
%
%    See also VERPINV, VERFULLCOLRANK, VERIFYLSS.

%    Copyright 2008 Jiri Rohn
%
%    Based on the description of the set X of all least squares solutions as
%        X = { pinv(A)*b+(eye(n,n)-pinv(A)*A)*y | y in R^n }
%    (see e.g. Stewart and Sun), with a verified pseudoinverse being computed
%    by VERPINV.
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
%    2008-03-20   first version
%    2008-03-31   version for posting
%    2008-04-05   output variable E added
%
gr=getround;
setround(0);
b=b(:); 
[m,n]=size(A);  
x=repmat(infsup(NaN,NaN),n,1);
B=repmat(infsup(NaN,NaN),n,n);
E.error='verlsq: none';
E.where='NaN';
E.value='NaN';
% data check
if nargin~=2||nargout>3||m~=length(b)||~isreal(A)||~isreal(b)||isintval(A)||isintval(b)
    E.error='verlsq: wrong data';
    setround(gr); return
end
% pseudoinverse
[X,Everpinv]=verpinv(A); % verified pseudoinverse; intval quantity
if isnan(X.inf(1,1)) % pseudoinverse not computed
    E=Everpinv;
    setround(gr); return
end
% computation of x
x=X*b; 
% computation of B
I=eye(n,n);
if verfullcolrank(A)~=1 % full column rank of A not verified
    B=I-X*A; 
else % full column rank of A verified, lsq solution unique
    B=repmat(infsup(0,0),n,n); % B set to zero interval matrix
end
setround(gr);
