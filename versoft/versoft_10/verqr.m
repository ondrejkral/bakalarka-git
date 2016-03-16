function [Q,R,E]=verqr(A)
%    VERQR      Verified QR decomposition of an interval (or real) matrix. 
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For an interval (or real) mxn matrix A, 
%        [Q,R,E]=verqr(A)
%    computes an interval matrix Q and an upper triangular interval matrix
%    R having the following property: for each Ao in A there holds 
%        Ao=Qo*Ro     
%    (in exact arithmetic) for some Qo in Q with orthogonal columns and for
%    some Ro in R with nonnegative diagonal entries (the QR decomposition
%    of Ao). Here, 
%        if m>=n, then Q is m-by-n and R is n-by-n,
%        if  m<n, then Q is m-by-m and R is m-by-n. 
%
%    If no verified result is found, then Q, R consist of NaN's (R upper
%    triangular). 
%
%    The structured array E explains reasons for NaN output. It has three
%    fields: E.error, E.where, E.value. 
%
%    WARNING. Output data widths may grow rapidly with increasing dimensions.
%
%    See also QR, CHOL.

%    Copyright 2008 Jiri Rohn.
%
%    Computed as A'*A=L*L' (Cholesky decomposition), R=L', Q=A*inv(R) for m>=n,
%    and using Householder reflectors for m<n.  
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
%    2008-03-12   version for posting
%    2008-03-28   error output E added
%    2008-03-29   triu(R) used
%    2008-03-30   B=hull(B,B') added, verhouseqr merged
%    2008-03-31   version for posting
%    2008-04-16   warning added
%    2008-04-17   default Q, R added (missing)
%
gr=getround;
setround(0);
[m,n]=size(A);
if m>=n
    Q=repmat(infsup(NaN,NaN),m,n);
    R=repmat(infsup(NaN,NaN),n,n);
    R=triu(R);
else % m<n
    Q=repmat(intval(NaN),m,m);
    R=repmat(intval(NaN),m,n);
    R=triu(R);
end
E.error='verqr: none';
E.where='NaN';
E.value='NaN';
if ~(nargin==1&&nargout<=3&&isreal(A))
    E.error='verqr: wrong data';
    setround(gr); return
end
if ~isintval(A)
    A=infsup(A,A);
end
if m>=n % via Cholesky (calls verchol)
    Q=repmat(infsup(NaN,NaN),m,n);
    R=repmat(infsup(NaN,NaN),n,n);
    R=triu(R); % default upper triangular matrix of NaN's
    B=A'*A;
    B=hull(B,B'); % to compensate errors created by B=A'*A (otherwise may be unsymmetric) 
    [ch,L,Everchol]=verchol(B);
    if ~isequal(ch,1) % Cholesky decomposition not found
        E=Everchol;
        setround(gr); return
    end
    R=L'; % R is nxn upper triangular with positive diagonal
    Ri=inv(R);
    if isnan(Ri.inf(1,1)) % inverse not computed
        E.error='verqr: inverse of R not computed';
        setround(gr); return
    end
    Q=A*Ri; % Q is mxn
    if ~(all(all(in(infsup(eye(n,n),eye(n,n)),Q'*Q)))&&all(diag(R.inf)>0)) % to be sure
        E.error='verqr: Q not orthogonal, or diagonal of R not positive';
        setround(gr); return
    end
else % m<n % via Householder (calls verhouseqr)
    [Q,R,Everhouseqr]=verhouseqr(A);
    if isnan(Q.inf(1,1)) % not computed
        E=Everhouseqr;
    end
end
setround(gr);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions verhouseqr, house
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [Q,R,E]=verhouseqr(A)
% A=Q*R, where Q is orthogonal and R upper triangular
% with nonnegative diagonal entries
if ~isintval(A)
    A=infsup(A,A);
end
E.error='verhouseqr: none';
E.where='NaN';
E.value='NaN';
[m,n]=size(A);
Q=eye(m,m); Q=intval(Q);
for j=1:min(m,n)
    a=A(j:m,j);
    [x,Ehouse]=house(a); % Householder vector
    if isnan(x.inf(1))
        Q=repmat(intval(NaN),m,m);
        R=repmat(intval(NaN),m,n);
        R=triu(R);
        E=Ehouse;
        E.where=['j = ',int2str(j)];
        return
    end
    A(j:m,:)=A(j:m,:)-2*x*(x'*A(j:m,:)); % main formulae % premultiplying by Householder matrix
    Q(j:m,:)=Q(j:m,:)-2*x*(x'*Q(j:m,:));
end
if ~(all(all(in(infsup(eye(m,m),eye(m,m)),Q'*Q)))) % to be sure
    E.error='verhouseqr: Q not orthogonal';
    setround(gr); return
end
Q=Q';
R=A;
R=triu(R); % zeroing subdiagonal entries
for j=1:min(m,n) % making diagonal of R nonnegative (wherever possible)
    if R.sup(j,j)<=0
        R(j,:)=-R(j,:);
        Q(:,j)=-Q(:,j);
    end
end
%
function [x,E]=house(a)
% interval Householder vector
% a, x interval vectors
% (I-2*x*x')*a is a multiple of I(:,1) (possibly negative), x normalized in 2-norm
E.error='house: none';
E.where='NaN';
E.value='NaN';
x=a(:); n=length(x);
if isequal(x(2:n),intval(zeros(n-1,1))) % a has the required form   
    x=intval(zeros(n,1)); return
end
if x.mid(1)>=0 % preventing cancellation
    x(1)=x(1)+norm(x,2);
else
    x(1)=x(1)-norm(x,2);
end
nx=norm(x,2);
if isequal(nx.inf,0)
    x=repmat(intval(NaN),n,1);
    E.error='house: normalizing not possible';
    E.value=['nx = [',num2str(nx.inf),', ',num2str(nx.sup),']'];
    return
end
x=x/nx; % normalizing

