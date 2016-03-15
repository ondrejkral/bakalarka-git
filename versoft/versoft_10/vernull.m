function [N,E]=vernull(A)
%    VERNULL        Verified null space of a rectangular real matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a rectangular m-by-n real matrix A,
%        [N,E]=vernull(A)
%    computes an n-by-n interval matrix N verified to contain a real matrix
%    No such that the null space 
%        Null(A) = { x | A*x=0 }
%    is described by
%        Null(A) = { No*y | y in R^n }
%    (i.e., the columns of No span Null(A)). If no verified result is
%    found, then N consists of NaN's.  
%           
%    The structure E explains reasons for NaN output. 
%
%    COMMENT. The size of No is always n-by-n, although it should be
%    n-by-(n-r), where r is the rank of A. Hence, the columns of No may be
%    linearly dependent. Unfortunately, I have not found a better way so
%    far.
%    
%    See also VERBASIS, VERRANK, VERLSQ.

%    Copyright 2008 Jiri Rohn.
%
%    Verified null space N computed by VERLSQ as N=eye(size(A,2))-pinv(A)*A.
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
[m,n]=size(A); 
N=repmat(intval(NaN),n,n);
E.error='vernull: none';
E.where='NaN';
E.value='NaN';
if (nargin~=1)||(nargout>2)||~isreal(A)||isintval(A)
    E.error='vernull: wrong data';
    setround(gr); return                                   
end
if issparse(A)
    A=full(A); % sparse not implemented 
end
b=ones(m,1);
[x,N1,Everlsq]=verlsq(A,b);
if ~isnan(N1.inf(1,1)) % computed
    N=N1;
else % not computed
    E=Everlsq; % verlsq error message
end
setround(gr);
