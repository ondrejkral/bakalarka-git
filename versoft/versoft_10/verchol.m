function [ch,L,E]=verchol(A) 
%    VERCHOL     Verified Cholesky decomposition of a symmetric interval (or real) matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a symmetric interval (or real) matrix A,
%        [ch,L,E]=verchol(A)
%    computes a verified Cholesky factor L of A, or yields no verified result: 
%
%    ch= 1           L is a lower triangular interval matrix with positive
%                    diagonal entries verified for each symmetric Ao in A
%                    to contain a matrix Lo satisfying Ao=Lo*Lo' (Cholesky
%                    decomposition of Ao); this also proves that each
%                    symmetric Ao in A is positive definite,
%    ch= 0           each symmetric Ao in A verified not to be positive definite
%                    (L is a lower triangular matrix of NaN's), 
%    ch=-1           no verified result (L is a lower triangular matrix of NaN's).
%
%    The structured array E explains reasons for NaN output. It has three
%    fields: E.error, E.where, E.value. 
%
%    WARNING. Output data widths may grow rapidly with increasing dimensions.
%
%    See also CHOL.

%    Copyright 2008 Jiri Rohn.
%
%    Classical Cholesky decomposition performed in interval arithmetic. 
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
%    2008-02-02   first version
%    2008-02-14   version for posting
%    2008-03-05   formula for L(i,k) changed
%    2008-03-11   added as a subroutine to verqr
%    2008-03-28   error output E added
%    2008-03-29   main formulae vectorized, tril(L) used
%    2008-03-30   L eliminated, computation done in frame of A; 
%                 made independent file again
%    2008-03-31   version for posting
%    2008-04-16   warning added
%    2008-04-17   default L added (missing)
%
gr=getround;
setround(0);
[m,n]=size(A);
ch=-1;
L=repmat(infsup(NaN,NaN),n,n); L=tril(L);
E.error='verchol: none';
E.where='NaN';
E.value='NaN';
if ~(nargin==1&&nargout<=3&&m==n&&isreal(A)&&isequal(A,A')) % wrong data
    E.error='verchol: wrong data';
    setround(gr); return
end
if ~isintval(A)
    A=infsup(A,A);
end
% columnwise computation of L % done in frame of A, 2008-03-30
for k=1:n
    el=(A(k,1:k-1))'; % column vector % enables vectorized computation
    alpha=A(k,k)-el'*el; % first main formula (diagonal entry)
    if alpha.inf<=0 
        if alpha.sup<=0
            ch=0; 
            L=repmat(infsup(NaN,NaN),n,n); L=tril(L); % each symmetric Ao in A verified not to be PD
            setround(gr); return
        else % alpha.inf<=0, alpha.sup>0
            E.error='verchol: pivot not positive';
            E.where=['k = ',int2str(k)];
            E.value=['alpha = [',num2str(alpha.inf),', ',num2str(alpha.sup),']'];
            ch=-1; 
            L=repmat(infsup(NaN,NaN),n,n); L=tril(L); % no verified result
            setround(gr); return
        end
    end
    % alpha.inf>0
    A(k,k)=sqrt(alpha);
    if A(k,k).inf<=0 % to be sure
        E.error='verchol: square root not positive';
        E.where=['k = ',int2str(k)];
        E.value=['L(k,k) = [',num2str(A(k,k).inf),', ',num2str(A(k,k).sup),']'];
        ch=-1; L=L0; % no verified result
        setround(gr); return
    end
    A(k+1:n,k)=A(k+1:n,k)-A(k+1:n,1:k-1)*el; % second main formula (subdiagonal entries) % changed 2008-03-29 to vectorized version
    A(k+1:n,k)=A(k+1:n,k)/A(k,k);
end
ch=1; % verified Cholesky decomposition found
L=tril(A); % lower triangular part extracted
setround(gr);
