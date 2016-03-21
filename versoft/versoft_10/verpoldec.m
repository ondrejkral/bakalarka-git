function [P,Q,E]=verpoldec(A)
%    VERPOLDEC        Verified polar decomposition of a square real matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a square real matrix A,
%        [P,Q,E]=verpoldec(A)
%    computes square interval matrices P, Q of the same size that are
%    verified to contain a positive semidefinite real matrix Po and an
%    orthogonal real matrix Qo satisfying
%        A=Po*Qo
%    in exact arithmetic (polar decomposition; Autonne 1902). If A is
%    nonsingular, then both Po, Qo are unique. If no verified result is
%    found, then P, Q consist of NaN's.  
%           
%    The structure E explains reasons for NaN output.
%   
%    See also VERTHINSVD.

%    Copyright 2008 Jiri Rohn.
%
%    Computed from thin SVD of A using VERTHINSVD.
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
P=repmat(intval(NaN),m,n);
Q=repmat(intval(NaN),m,n);
E.error='verpoldec: none';
E.where='NaN';
E.value='NaN';
% data check
if ~(nargin==1&&nargout<=3&&m==n&&isreal(A)&&~isintval(A)) % wrong data
    E.error='verpoldec: wrong data';
    setround(gr); return
end
% thin svd
[U,S,V,Everthinsvd]=verthinsvd(A);
if isnan(U.inf(1,1)) % thin SVD not computed
    E=Everthinsvd;
    setround(gr); return
end
% thin SVD computed: A=U*S*V'; U, S, V: nxn
P=U*S*U'; % positive semidefinite
Q=U*V'; % orthogonal
% P*Q=U*S*U'*U*V'=U*S*V'=A
setround(gr);
