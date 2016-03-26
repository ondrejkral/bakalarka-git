function [X,E]=vermatreqn(A,B,C,D,F)
%    VERMATREQN       Verified solution of the matrix equation A*X*B+C*X*D=F
%                     (in particular, of the Sylvester or Lyapunov equation).   
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For interval (or real) m-by-m matrices A, C, p-by-p matrices B, D, and
%    an m-by-p matrix F,
%        [X,E]=vermatreqn(A,B,C,D,F)
%    returns a verified m-by-p solution X of the matrix equation
%        A*X*B+C*X*D=F,     
%    or yields no verified result.
%
%    Possible outcomes:
%
%    ~isnan(X.inf(1,1)) :  it is verified that for each Ao in A, Bo in B,
%                          Co in C, Do in D, and Fo in F, the equation
%                             Ao*Xo*Bo+Co*Xo*Do=Fo
%                          possesses a unique matrix solution Xo which is
%                          verified to be contained in X,
%     isnan(X.inf(1,1)) :  no verified result (the interval matrix X
%                          consists of NaN's).  
%
%    The structure E explains reasons for NaN output.
%
%    For a verified solution of the SYLVESTER equation
%        A*X+X*D=F
%    (A, D square of possibly different sizes), use 
%        [X,E]=vermatreqn(A,eye(size(D)),eye(size(A)),D,F).
%
%    For a verified solution of the LYAPUNOV equation
%        A*X+X*A'=F 
%    (A, F square of the same size), use 
%        [X,E]=vermatreqn(A,eye(size(A')),eye(size(A)),A',F).
%
%    See also VERIFYLSS.

%    Copyright 2008 Jiri Rohn
%
%    Based on solving the square system of linear equations
%        (kron(B',A)+kron(D',C))*vec(X)=vec(F),
%    where vec(X)=X(:) and vec(F)=F(:); see R. A. Horn and C. R. Johnson,
%    Topics in Matrix Analysis, Cambridge University Press, Cambridge 1991,
%    Lemma 4.3.1.  
%
%    This work was done during author's employment at the Anglo-American
%    University in Prague, Czech Republic.
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
%    2008-05-16   first version 
%    2008-05-27   version for posting
%
gr=getround;
setround(0);
[m,n]=size(A); [p,q]=size(B);
% defaults
X=repmat(infsup(NaN,NaN),n,p);
E.error='vermatreqn: none';
E.where='NaN';
E.value='NaN';
% data check
if ~(nargin==5&&nargout<=2&&m==n&&p==q&&all(size(A)==size(C))&&all(size(B)==size(D))&&all(size(F)==[m,q]))
    E.error='vermatreqn: wrong data';
    setround(gr); return
end
% equation:  A*X*B+C*X*D=F % A, C: mxm, B, D: pxp
% solved as: (kron(B',A)+kron(D',C))*vec(X)=vec(F) % see Horn and Johnson, Topics ..., Lemma 4.3.1
% the system matrix is mpxmp
AK=kronmod(B',A)+kronmod(D',C);
if issparse(AK)
    AK=full(AK); % sparse matrices not implemented in verifylss
end
bK=F(:);
xK=verifylss(AK,bK); % matrix in vector form
X=reshape(xK,n,p);   % reshaped back
if isnan(X.inf(1,1)) % no result
    E.error='vermatreqn: verifylss failed';
end
setround(gr);
%
%%%%%%%%%%%%%%%%%%%%%
% Subfunction kronmod
%%%%%%%%%%%%%%%%%%%%%
%
function K=kronmod(A,B)
% Kronecker product of matrices of type intval
if ~isintval(A)
    A=infsup(A,A);
end
if ~isintval(B)
    B=infsup(B,B);
end
% the following five lines (construction of the Kronecker product) 
% copied from KRON.M, Copyright 1984-2004 The MathWorks, Inc. 
% (but the last product is performed in interval arithmetic)
[ma,na]=size(A);
[mb,nb]=size(B);
[ia,ib]=meshgrid(1:ma,1:mb);
[ja,jb]=meshgrid(1:na,1:nb);
K=A(ia,ja).*B(ib,jb);

