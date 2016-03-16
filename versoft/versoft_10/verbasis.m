function [B,K,E]=verbasis(A)
%    VERBASIS        Verified basis of the range space of a rectangular real matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a rectangular real matrix A,
%        [B,K,E]=verbasis(A)
%    computes an index set K such that B=A(:,K) is a verified basis of the
%    range space of A. (Thus, B is real, not interval matrix, and
%    r=length(K) is the verified rank of A.) If no verified result is
%    found, then B, K consist of NaN's. 
%           
%    The structure E explains reasons for NaN output.  
%    
%    See also RREF, VERFULLCOLRANK, RANK, VERRANK.

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
%    2008-04-12   first version 
%    2008-04-18   version for posting
%
gr=getround;
setround(0);
[m,n]=size(A); 
B=repmat(NaN,m,n);
K=repmat(NaN,1,n);
E.error='verbasis: none';
E.where='NaN';
E.value='NaN';
if (nargin~=1)||(nargout>3)||~isreal(A)||isintval(A)
    E.error='verbasis: wrong data';
    setround(gr); return                                   
end
if issparse(A)
    A=full(A); % sparse not implemented 
end
[AR,K]=rref(A); % K index set of the expected basis (later B)
if isempty(K) % no basis index set found
    K=repmat(NaN,1,n);
    E.error='verbasis: no basis index set found';
    setround(gr); return
end
% K nonempty
B=A(:,K); % expected basis
fcr=verfullcolrank(B);
if fcr~=1 % columns of B not verified linearly independent
    B=repmat(NaN,m,n);
    K=repmat(NaN,1,n);
    E.error='verbasis: columns of expected basis B not verified linearly independent';
    setround(gr); return
end
% columns of B verified linearly independent
for j=1:n
    if ~any(K==j)
        fcr=verfullcolrank([B A(:,j)]);
        if fcr~=0 % columns of [B A(:,j)] not verified linearly dependent
            B=repmat(NaN,m,n);
            K=repmat(NaN,1,n);
            E.error='verbasis: j-th column of A not verified to belong to the range space of B';
            E.where=['j = ',int2str(j)];
            setround(gr); return
        end
        % A(:,j) verified to belong to the range space of B
    end
end
% B verified to span the range space of A
setround(gr);















