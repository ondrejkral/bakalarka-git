function [y,E]=verorthproj(A,x)
%    VERORTHPROJ        Verified orthogonal projection on the range space of a rectangular real matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a rectangular m-by-n real matrix A and an m-by-1 real vector x,
%        [y,E]=verorthproj(A,x)
%    computes an interval vector y which is verified to contain the
%    orthogonal projection yo of the vector x on the range space R(A) of A
%    (i.e., yo belongs to R(A) and has minimal Euclidean distance from x).
%    If no verified result is found, then y consists of NaN's.  
%
%    The structure E explains reasons for NaN output. 
%
%    See also PINV, VERPINV.

%    Copyright 2008 Jiri Rohn.
%
%    Verified projection computed by VERPINV as y=A*pinv(A)*x.
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
%    2008-04-13   first version
%    2008-04-18   version for posting
%
gr=getround;
setround(0);
m=size(A,1); 
y=repmat(intval(NaN),m,1);
E.error='verorthproj: none';
E.where='NaN';
E.value='NaN';
if (nargin~=2)||(nargout>2)||~isreal(A)||isintval(A)||length(x)~=m
    E.error='verorthproj: wrong data';
    setround(gr); return                                   
end
if issparse(A)
    A=full(A); % sparse not implemented 
end
[X,Everpinv]=verpinv(A);
if ~isnan(X.inf(1,1)) % pseudoinverse computed
    if verfullcolrank(A')==1 % identical mapping
        y=infsup(x,x); 
    else % projection
        y=A*X*x; 
    end
else % not computed
    E=Everpinv; % verpinv error message
end
setround(gr);

