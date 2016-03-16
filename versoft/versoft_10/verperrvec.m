function [pv,rho,As]=verperrvec(A,x)
%    VERPERRVEC      Verified Perron vector of a positive interval matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a positive square interval matrix A and a positive vector x,
%    [pv,rho,As]=verperrvec(A,x)
%    verifies x/norm(x,1) to be the Perron vector of some matrix in A,
%    or not to be the Perron vector of any matrix in A, or yields no
%    verified result:  
%
%    pv= 1           x/norm(x,1) is verified to be the Perron vector of some matrix in A,
%                    rho is an interval number such that for each rho0 in rho, 
%                       A is verified to contain a matrix having spectral radius rho0 
%                       and Perron vector x/norm(x,1), 
%                    As is a very tight interval matrix verified to be a part of A and 
%                       to contain a matrix having spectral radius mid(rho)
%                       and Perron vector x/norm(x,1),  
%    pv= 0           x is verified not to be the Perron vector of any matrix in A,
%                       rho and As consist of NaN's, 
%    pv=-1           no verified result (data may be wrong).
%
%    COMMENT. The Perron vector of a positive real matrix Ao is the unique
%    positive eigenvector x pertaining to the positive eigenvalue
%    rho(Ao)=max(abs(eig(Ao))) and satisfying sum(x)=1. Since the equality
%    sum(x)=1 cannot be generally verified in floating point, we circumvent this
%    difficulty by verifying that the positive vector x, after being
%    normalized by x/norm(x,1) (in infinite precision), will become the
%    Perron vector of Ao. Alternatively, we can see it as verified that x is
%    a positive multiple of the Perron vector of Ao.
%
%    See also VEREIGVEC.

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
gr=getround;
setround(0);
% checking data
x=x(:);
[m,n]=size(A); p=length(x);
pv=-1; rho=infsup(NaN,NaN); As=repmat(infsup(NaN,NaN),m,n);
if m~=n||n~=p||~isreal(A)||~isintval(A)||~isreal(x)||isintval(x)||~all(all(A.inf>0))||~all(x>0) % error('wrong data')
   setround(gr); return
end
% main part: application of vereigvec
[pv,rho,As]=vereigvec(A,x);
setround(gr);
















