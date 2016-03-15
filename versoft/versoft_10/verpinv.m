function [X,E]=verpinv(A)
%    VERPINV      Verified pseudoinverse of a real matrix. 
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a rectangular real matrix A,
%        X=verpinv(A)
%    either computes a verified pseudoinverse X of A, or yields no verified
%    result in which case X consists of NaN's. 
%
%    The structured array E explains reasons for NaN output. It has three entries, 
%    each having three fields: E(j).error, E(j).where, E(j).value, j=1:3.  
%
%    See also INV, VERINVERSE, VERIFYLSS.

%    Copyright 2008 Jiri Rohn
%
%    Computed as intersection of enclosures obtained first by using explicit
%    formulae (in case of full rank), second by employing Greville's
%    algorithm; if not successful, third by using SVD. 
%
%    This work was partly supported by the Czech Republic National Research
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
%    History
%
%    2007-10-07   first version (all subroutines included)
%    2008-01-14   version for posting
%    2008-04-04   output variable E added
%    2008-04-05   verpinvviathinsvd added, former verpinv renamed as verpinvold
%
% set rounding to nearest
gr=getround;
setround(0);
[m,n]=size(A); 
X=repmat(infsup(NaN,NaN),n,m);
E.error='verpinv: none';
E.where='NaN';
E.value='NaN';
if (nargin~=1)||~isreal(A)||isintval(A)
    E.error='verpinv: wrong data';
    setround(gr); return                                   
end
if issparse(A)
    A=full(A); % sparse not implemented 
end
% two subroutines
[X1,Everpinvold]=verpinvold(A);
if isnan(X1.inf(1,1)) % this "if" added later; three enclosures too time-consuming
    [X2,Everpinvviathinsvd]=verpinvviathinsvd(A);
else % enclosure computed
    X2=repmat(infsup(NaN,NaN),n,m);
end
% output cases
if ~isnan(X1.inf(1,1))&&~isnan(X2.inf(1,1)) % does not occur with the above "if"; left for possible later changes
    X=intersect(X1,X2);
    setround(gr); return
else
    if ~isnan(X1.inf(1,1))
        X=X1;
        setround(gr); return
    end
    if ~isnan(X2.inf(1,1))
        X=X2;
        setround(gr); return
    end
    % no enclosure
    clear E; % now structured array
    E=Everpinvold; % E(1), E(2) contained
    E(3)=Everpinvviathinsvd; 
end
setround(gr);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions: verpinvold, verpinvsimple, verpinvbasic, vergreville, vergrevillebasic, verpinvviathinsvd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [X,E]=verpinvold(A) % former verpinv (up to 2008-04-05)
gr=getround;
setround(0);
[m,n]=size(A); 
X=repmat(infsup(NaN,NaN),n,m);
E.error='verpinvold: none';
E.where='NaN';
E.value='NaN';
if (nargin~=1)||~isreal(A)||isintval(A)
    E.error='verpinvold: wrong data';
    setround(gr); return                                   
end
if issparse(A)
    A=full(A); % sparse not implemented 
end
% two subroutines
[X1,Everpinvsimple]=verpinvsimple(A);
[X2,Evergreville]=vergreville(A);
% output cases
if ~isnan(X1.inf(1,1))&&~isnan(X2.inf(1,1))
    X=intersect(X1,X2);
    setround(gr); return
else
    if ~isnan(X1.inf(1,1))
        X=X1;
        setround(gr); return
    end
    if ~isnan(X2.inf(1,1))
        X=X2;
        setround(gr); return
    end
    clear E; % now structured array
    E(1)=Everpinvsimple;
    E(2)=Evergreville; 
end
setround(gr);
%
function [X,E]=verpinvsimple(A)
% computes verified pseudoinverse under full column or row rank
[m,n]=size(A);
E.error='verpinvsimple: none';
E.where='NaN';
E.value='NaN';
if m>=n
    [X,Everpinvbasic]=verpinvbasic(A);
    E=Everpinvbasic;
else
    [X,Everpinvbasic]=(verpinvbasic(A'));
    X=X';
    E=Everpinvbasic;
end
%
function [X,E]=verpinvbasic(A)
% computes verified pseudoinverse under full column rank by explicit formulae
[m,n]=size(A);
A=infsup(A,A); % intval quantity
X=repmat(infsup(NaN,NaN),n,m);
E.error='verpinvbasic: none';
E.where='NaN';
E.value='NaN';
if m<n
    E.error='verpinvbasic: wrong data: m<n';
    return
end
% two ways
X1=inv(A'*A)*A';
I=eye(m,m); I=infsup(I,I); % mxm
O=zeros(n,n); O=infsup(O,O); % nxn
Onm=zeros(n,m); Onm=infsup(Onm,Onm); % nxm
A1=[A -I; O A']; % (m+n)x(n+m)
B1=[I; Onm]; % (m+n)xm
X2=verifylss(A1,B1); % system: A*X-Y=I, A'*Y=0
X2=X2(1:n,:); % X part
% output cases
if ~isnan(X1.inf(1,1))&&~isnan(X2.inf(1,1))
    X=intersect(X1,X2);
    return
else
    if ~isnan(X1.inf(1,1))
        X=X1;
        return
    end
    if ~isnan(X2.inf(1,1))
        X=X2;
        return
    end
    E.error='verpinvbasic: none of the two applications of verifylss successful';
end
%
function [X,E]=vergreville(A)
% computes verified pseudoinverse X of A by Greville's algorithm 
% applied to both A and A'
[m,n]=size(A); 
X=repmat(infsup(NaN,NaN),n,m);
E.error='vergreville: none';
E.where='NaN';
E.value='NaN';
[X1,Evergrevillebasic]=vergrevillebasic(A);
X2=(vergrevillebasic(A'))';
% output cases
if ~isnan(X1.inf(1,1))&&~isnan(X2.inf(1,1))
    X=intersect(X1,X2);
    return
else
    if ~isnan(X1.inf(1,1))
        X=X1;
        return
    end
    if ~isnan(X2.inf(1,1))
        X=X2;
        return
    end
    E=Evergrevillebasic; % only one error can be outputed
end
%
function [X,E]=vergrevillebasic(A)
% computes verified pseudoinverse X of A by Greville's algorithm
[m,n]=size(A); 
E.error='vergrevillebasic: none';
E.where='NaN';
E.value='NaN';
A=infsup(A,A);
d=A(:,1);
if isequal(d,repmat(infsup(0,0),m,1))
    X=d';  
else
    scal=d'*d;
    if scal.inf>0
        X=d'/scal;
    else
        X=repmat(infsup(NaN,NaN),n,m);
        E.error='vergrevillebasic: division by zero';
        E.where=['j = ',int2str(1)];
        E.value=['scal = [',num2str(scal.inf),', ',num2str(scal.sup),']'];
        return
    end
end
for j=2:n
    d=X*A(:,j);
    c=A(:,j)-A(:,1:(j-1))*d;
    if ~(any(c.sup<0)||any(c.inf>0))
        X=repmat(infsup(NaN,NaN),n,m);
        E.error='vergrevillebasic: neither c==0 nor c~=0 could be verified';
        E.where=['j = ',int2str(j)];
        E.value=['c = [',num2str(c.inf'),', ',num2str(c.sup'),']'];
        return
    end
    scal=c'*c;
    if scal.inf>0
        bt=c'/(c'*c);
    else
        X=repmat(infsup(NaN,NaN),n,m);
        E.error='vergrevillebasic: division by zero';
        E.where=['j = ',int2str(j)];
        E.value=['scal = [',num2str(scal.inf),', ',num2str(scal.sup),']'];
        return
    end
    X=[X-d*bt; bt]; % Greville's formula
end
%
function [X,E]=verpinvviathinsvd(A)
% computes verified pseudoinverse X of A via thin SVD
[m,n]=size(A); 
X=repmat(infsup(NaN,NaN),n,m);
E.error='verpinvviathinsvd: none';
E.where='NaN';
E.value='NaN';
[U,S,V,Everpinvviathinsvd]=verthinsvd(A); % A=U*S*V'
if isnan(U.inf(1,1)) % thin SVD not computed
    E=Everpinvviathinsvd;
    return
end
% thin SVD computed
sigma=diag(S); % interval vector of singular values % S: nxn
r1=length(find(sigma.inf>0)); % number of positive lower bounds
r2=length(find(sigma.sup>0)); % number of positive upper bounds
if r1~=r2
    r=infsup(r1,r2);
    E.error='verpinvviathinsvd: verified rank could not be established';
    E.value=['r = [',int2str(r.inf),', ',int2str(r.sup),']'];
end
r=r1; % verified rank
if r>0 % S~=0
    sigma=sigma(1:r); % positive singular values
    sigma=1./sigma; % inverse
    S(1:r,1:r)=diag(sigma);
    X=V*S*U'; % main formula for pseudoinverse
else % S=0, A=0
    X=repmat(infsup(0,0),n,m);
end
