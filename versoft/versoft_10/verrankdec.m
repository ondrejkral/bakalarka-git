function [B,C,E]=verrankdec(A)
%    VERRANKDEC        Verified rank decomposition of a rectangular real matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a rectangular m-by-n real matrix A,
%        [B,C,E]=verrankdec(A)
%    computes a real m-by-r matrix B and an interval r-by-n matrix C
%    verified to contain a real matrix Co such that
%        A=B*Co
%    holds in exact arithmetic. Moreover, r is the verified rank of A, B has
%    linearly independent columns and Co has linearly independent rows (a
%    rank decomposition of A). If no verified result is found, then B, C
%    consist of NaN's.  
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
%    2008-04-12   first version 
%    2008-04-18   version for posting
%
gr=getround;
setround(0);
[m,n]=size(A); 
B=repmat(NaN,m,n);
C=repmat(NaN,n,n);
E.error='verrankdec: none';
E.where='NaN';
E.value='NaN';
if (nargin~=1)||(nargout>3)||~isreal(A)||isintval(A)
    E.error='verrankdec: wrong data';
    setround(gr); return                                   
end
if issparse(A)
    A=full(A); % sparse not implemented 
end
[AR,K]=rref(A); % K index set of the expected basis (later B)
if isempty(K) % no basis index set found
    E.error='verrankdec: no basis index set found';
    setround(gr); return
end
% K nonempty
B=A(:,K); % expected basis
fcr=verfullcolrank(B);
if fcr~=1 % columns of B not verified linearly independent
    B=repmat(NaN,m,n);
    E.error='verrankdec: columns of expected basis B not verified linearly independent';
    setround(gr); return
end
% columns of B verified linearly independent
[Bp,Everpinv]=verpinv(B);
if isnan(Bp.inf(1,1)) % pseudoinverse of B not computed
    B=repmat(NaN,m,n);
    E=Everpinv;
    setround(gr); return
end
% pseudoinverse of B computed (intval quantity)
r=length(K); 
C=repmat(intval(0),r,n);
Ir=intval(eye(r,r));
for j=1:n
    if ~any(K==j) % j not in K
        fcr=verfullcolrank([B A(:,j)]);
        if fcr~=0 % columns of [B A(:,j)] not verified linearly dependent
            B=repmat(NaN,m,n);
            C=repmat(NaN,n,n);
            E.error='verrankdec: j-th column of A not verified to belong to the range space of B';
            E.where=['j = ',int2str(j)];
            setround(gr); return
        end
        % A(:,j) verified to belong to the range space of B
        C(:,j)=Bp*A(:,j); % intval quantity
    else % j in K
        i=find(K==j); % j at the i-th place of K
        C(:,j)=Ir(:,i);
    end
end
% first enclosure C found
% as a by-product, r is the verified rank of A
% second enclosure
B0=intval(B);
BTBi=inv(B0'*B0);
if isnan(BTBi.inf(1,1)) % inverse not computed
    E.where='only first enclosure outputed';
    setround(gr); return % first enclosure outputed
end
% inverse computed
C1=(BTBi*B')*A; % A=B*C implies B'*A=B'*B*C implies C=inv(B'*B)*B'*A
C2=intersect(C,C1); % intersection of the two enclosures
if any(any(isnan(C2))) % to be sure
    E.where='only first enclosure outputed';
    setround(gr); return % first enclosure outputed
else
    C=C2; % intersection outputed
end
% verified A=B*Co, Co in C, B having linearly independent columns
setround(gr);
