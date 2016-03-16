function [B,As]=verinverse(A)
%    VERINVERSE     Verified inverse interval matrix.
%
%    For a square interval matrix A,
%    [B,As]=verinverse(A)
%    either computes a verified interval inverse B of A, 
%    or finds a tight ("almost thin") interval matrix As which is 
%    a part of A and is verified to contain a real singular matrix, 
%    or yields no verified result.
%
%    At least one of the interval matrices B, As always consists of NaN's
%    (as the two options exclude each other).
%
%    See also INV, VERINTERVALHULL.

%    Copyright 2007 Jiri Rohn
%
%    Based on the algorithm "inverse" described in
%    J. Rohn, A Handbook of Results on Interval Linear Problems,
%    posted at http://www.cs.cas.cz/~rohn
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
% set rounding to nearest
gr=getround;
setround(0);
% setting default output
[m,n]=size(A);
B=repmat(infsup(NaN,NaN),m,n); As=repmat(infsup(NaN,NaN),m,n); 
% checking compatibility of data
if (m~=n)||~isreal(A)
    setround(gr); return                                   % matrix not square or not real
end
% handling thin data
if ~isintval(A)                                            % allows for real input
    B=inv(infsup(A,A));                                    % inverse computed directly
    setround(gr); return
end
if issparse(A)
    A=full(A);                                             % sparse matrices not implemented yet
end
I=eye(n,n); I=infsup(I,I);
B=I;
% columnwise computation of the inverse
for j=1:n
    [xx,AAs]=intervalhull(A,I(:,j));                               
    if ~isempty(xx)                       
        B(:,j)=xx;                                         % verified j-th column computed
    end     
    if  ~isempty(AAs) 
        B=repmat(infsup(NaN,NaN),n,n);
        As=AAs;                                            % verified singular
        setround(gr); return           
    end
    if  isempty(xx)&&isempty(AAs)                              
        B=repmat(infsup(NaN,NaN),n,n);                     % no verified result
        setround(gr); return  
    end 
end
setround(gr);













