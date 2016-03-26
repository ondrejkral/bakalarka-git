function [Q,E]=verorth(A)
%    VERORTH        Verified orthonormal basis of the range space of a rectangular real matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a rectangular real matrix A,
%        [Q,E]=verorth(A)
%    computes an interval matrix Q which is verified to contain a real
%    matrix Qo whose columns form an orthonormal basis for the range space
%    of A (thus, size(Q,2) is the verified rank of A). If no verified
%    result is found, then Q consists of NaN's.  
%           
%    The structure E explains reasons for NaN output.  
%    
%    See also VERBASIS, VERTHINSVD, ORTH.

%    Copyright 2008 Jiri Rohn.
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
%    2008-04-18   version for posting
%
gr=getround;
setround(0);
[m,n]=size(A); 
% defaults
Q=repmat(NaN,m,n);
E.error='verorth: none';
E.where='NaN';
E.value='NaN';
% data check
if (nargin~=1)||(nargout>2)||~isreal(A)||isintval(A)
    E.error='verorth: wrong data';
    setround(gr); return                                   
end
if issparse(A)
    A=full(A); % sparse not implemented 
end
% basis
[B,K,Everbasis]=verbasis(A);
if isnan(B(1,1)) % no verified basis found
   E=Everbasis;
   setround(gr); return
end
% verified basis found; B is m-by-r, m>=r, r rank of A
[U,S,V,Everthinsvd]=verthinsvd(B);
if isnan(U.inf(1,1)) % verified thin svd not computed
    E=Everthinsvd;
    setround(gr); return
end
% verified thin svd computed
s=diag(S);
if ~all(s.inf>0) % to be sure
    E.error='verorth: U not verified to be a basis';
    setround(gr); return
end
% S*V' nonsingular, U proved to be a basis for R(A)
Q=U; % verified orthonormal basis
setround(gr);
