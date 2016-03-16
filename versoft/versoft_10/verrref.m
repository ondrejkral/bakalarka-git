function [AR,E]=verrref(A)
%    VERRREF        Verified reduced row-echelon form (RREF) of a rectangular real matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a rectangular real matrix A,
%        [AR,E]=verrref(A)
%    computes a verified reduced row-echelon form (RREF) AR of A. In
%    particular, the number of nonzero rows of AR yields the verified rank
%    of A. The result is not computed by Gaussian elimination, which
%    contributes to its accuracy. If no verified result is found, then AR
%    consists of NaN's.  
%           
%    The structure E explains reasons for NaN output.  
% 
%    See also RREF, VERFULLCOLRANK, RANK, VERRANK, VERPINV.

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
%    2008-04-13   first version 
%    2008-04-19   version for posting
%
gr=getround;
setround(0);
[m,n]=size(A); 
AR=repmat(intval(NaN),m,n);
E.error='verrref: none';
E.where='NaN';
E.value='NaN';
if (nargin~=1)||(nargout>2)||~isreal(A)||isintval(A)
    E.error='verrref: wrong data';
    setround(gr); return                                   
end
if issparse(A)
    A=full(A); % sparse not implemented 
end
[AR,K]=rref(A); % K index set of the expected basis (later B)
K=sort(K); % increasing order, to be sure
if isempty(K) % no basis index set found
    AR=repmat(intval(NaN),m,n);
    E.error='verrref: no basis index set found';
    setround(gr); return
end
% K nonempty
B=A(:,K); % expected basis
fcr=verfullcolrank(B);
if fcr~=1 % columns of B not verified linearly independent
    AR=repmat(intval(NaN),m,n);
    E.error='verrref: columns of expected basis B not verified linearly independent';
    setround(gr); return
end
% columns of B verified linearly independent
r=length(K); 
AR=repmat(intval(0),m,n);
Ir=intval(eye(r,r));
for j=1:n
    if ~any(K==j) % j not in K
        B0=A(:,K(find(K<j))); % repaired 2008-04-19
        fcr=verfullcolrank([B0 A(:,j)]);
        if fcr~=0 % columns of [B0 A(:,j)] not verified linearly dependent
            AR=repmat(intval(NaN),m,n);
            E.error='verrref: j-th column of A not verified to belong to the range space of B0';
            E.where=['j = ',int2str(j)];
            setround(gr); return
        end
        % A(:,j) verified to belong to the range space of B0
        B0p=verpinv(B0);
        if isnan(B0p.inf(1,1)) % pseudoinverse of B0 not computed
            AR=repmat(intval(NaN),m,n);
            E.error='verrref: pseudoinverse of B0 not computed';
            E.where=['j = ',int2str(j)];
            setround(gr); return
        end
        % pseudoinverse of B0 computed (intval quantity)
        C=B0p*A(:,j); % intval quantity
        AR(1:length(C),j)=C;
    else % j in K
        i=find(K==j); % j at the i-th place of K
        C=Ir(:,i);
        AR(1:length(C),j)=C;
    end
end
% first enclosure AR found
% as a by-product, r is the verified rank of A
% second enclosure
C=AR(1:r,:); % first r rows of AR (A=B*C forms rank decomposition now)
B=intval(B);
BTBi=inv(B'*B);
if isnan(BTBi.inf(1,1)) % inverse not computed
    E.where='only first enclosure outputed';
    setround(gr); return % first enclosure outputed
end
% inverse computed
C1=BTBi*B'*A; % A=B*C implies B'*A=B'*B*C implies C=inv(B'*B)*B'*A
C2=intersect(C,C1); % intersection of the two enclosures
if any(any(isnan(C2))) % to be sure
    E.where='only first enclosure outputed';
    setround(gr); return % first enclosure outputed
else
    AR(1:r,:)=C2; % intersection outputed
end
% AR verified RREF of A
setround(gr);
