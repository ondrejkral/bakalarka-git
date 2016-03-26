function [pm,E,C]=verpmat(A,t,d)
%    VERPMAT        Verified P-property of a real matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a square real matrix A,
%        [pm,E,C]=verpmat(A)
%    either verifies that A is a P-matrix, or verifies that it is not a
%    P-matrix (in which case it yields a certificate of the fact), or
%    fails (yields an interval of NaN's).  
%
%    Possible outputs:
%
%    pm= 1:    A is verified to be a P-matrix,
%    pm= 0:    A is verified not to be a P-matrix; in this case either 
%              C.ij=i and C.submat=A(i,i)<=0, or C.ij=[i j] and the principal 
%              2-by-2 submatrix C.submat=[A(i,i) A(i,j); A(j,i) A(j,j)] 
%              has a verified nonpositive determinant C.det (a certificate),
%    pm=-1:    no verified result.
%
%    If pm~=0, then the C.ij, C.submat, C.det consist of NaN's.
%
%    The structured array E explains reasons for pm==-1. It has three fields: 
%    E.error, E.where, E.value.
%
%    The command
%        [pm,E,C]=verpmat(A,t)
%    (i.e., with an additional input argument t) performs the same task, but
%    in case of t==1 it also produces screen output of the form 
%        Expected remaining time: ... sec.
%    This feature, however, may slow down the actual computation. Use t=0
%    if you wish to suppress this feature. 
%
%    The command
%        [pm,E,C]=verpmat(A,t,d)
%    (i.e., with yet another additional input argument d) performs the same task,
%    but in case of d==1 it also describes on the screen which part of the
%    program is currently being executed, as e.g.
%       verifying not-P-property ...
%          not verified
%       verifying positive definiteness ...
%          not verified
%       verifying nonsingularity of A-I, A+I ...
%          verified
%       verifying P-property ...
%          computing the inverse of A-I ...
%          verifying regularity via verregsing ...
%          verified
%    Use d=0 if you wish to suppress this feature.
%
%    Checking not-P-property (pm==0) is performed in O(n^2) time, where
%    n=size(A,1) (it is based on a sufficient, but not necessary condition); 
%    however, checking the P-property (pm==1) occasionally may last long
%    since the problem is co-NP-complete (Coxson, Math. Progg. 1994). 
%
%    See also VERREGSING, VERPOSDEF.

%    Copyright 2008 Jiri Rohn.
%
%    Based on the following theorem by S. M. Rump in "On P-Matrices", Linear
%    Algebra and its Applications 363(2003), 237–250: if 
%        det(I-A)*det(I+A)~=0, 
%    then A is a P-matrix if and only if the interval matrix 
%        infsup(inv(A-I)*(A+I)-I,inv(A-I)*(A+I)+I) 
%    is regular, where I is the identity matrix. Regularity is checked with
%    the help of VERREGSING.
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
%    2008-11-26   first version
%    2008-11-27   A=1.1*A discarded, renamed (AA=intval(A)); use of veregsing for checking eigenvalues; t, d, verpd added
%    2008-12-01   help expanded
%    2008-12-03   final version (html)
%
gr=getround;
setround(0);
[m,n]=size(A);
% defaults
pm=-1;
E.error='verpmat: none';
E.where='NaN';
E.value='NaN';
C.ij=[NaN NaN];
C.submat=[NaN NaN; NaN NaN];
C.det=infsup(NaN,NaN);
% data check
if ~(nargin<=3&&nargout<=3&&m==n&&isreal(A)&&~isintval(A)) % wrong data
    E.error='verpmat: matrix not square or not real or of type intval';
    setround(gr); return
end
if ~(nargin>=2&&t==1)
    t=0;
end
if ~(nargin==3&&isequal(d,1))
    d=0;
end
AA=infsup(A,A); % A real, AA=infsup(A,A) interval 
% not-P-property case (usually faster)
if d==1
    disp('verifying not-P-property ...')
end
for i=1:n
    if A(i,i)<=0 % 1x1 minor nonpositive
        pm=0;
        C.ij=i;
        C.submat=A(i,i);
        C.det=intval(A(i,i));
        if d==1
            disp('   verified')
        end
        setround(gr); return % A is verified not to be a P-matrix
    end
    for j=i+1:n % i<j
        a=AA(i,i)*AA(j,j)-AA(i,j)*AA(j,i); % in interval arithmetic
        if a.sup<=0 % 2x2 minor nonpositive
            pm=0;
            C.ij=[i j];
            C.submat=[A(i,i) A(i,j)    
                      A(j,i) A(j,j)]; 
            C.det=a;
            if d==1
                disp('   verified')
            end
            setround(gr); return % A is verified not to be a P-matrix
        end
    end
end
if d==1
    disp('   not verified')
end
% checking positive definiteness
if d==1
    disp('verifying positive definiteness ...')
end
if verpd(A)==1
    pm=1;
    if d==1
        disp('   verified')
    end
    setround(gr); return % with verified P-property
end
if d==1
    disp('   not verified')
end
% checking regularity of A-I, A+I
if d==1
    disp('verifying nonsingularity of A-I, A+I ...')
end
I=eye(n,n); I=infsup(I,I);
if verregsing(AA-I)~=1||verregsing(AA+I)~=1 % checks via Gerschgorin circles or verified eigenvalues not proved effective
    E.error='verpmat: A-I, A+I not proved nonsingular; condition not verified';
    if verregsing(AA-I)~=1
        E.where='A-I';
        E.value=AA-I;
    else % verregsing(AA+I)~=1
        E.where='A+I';
        E.value=AA+I;
    end
    if d==1
        disp('   not verified')
    end
    setround(gr); return
end
if d==1
    disp('   verified')
end
% condition det(I-A)*det(I+A)~=0 satisfied
% regularity/P-property case (usually slower)
if d==1
    disp('verifying P-property ...')
end
if d==1
    disp('   computing the inverse of A-I ...')
end
B=inv(AA-I); % B of type intval
if isnan(B.inf(1,1))
    E.error='verpmat: verified inverse of A-I not computed';
    if d==1
        disp('   not verified')
    end
    setround(gr); return
end
% B computed
B=B*(AA+I); % B=inv(A-I)*(A+I)
Blow=B-I;
Bupp=B+I;
B=infsup(Blow.inf,Bupp.sup); % contains infsup(inv(A-I)*(A+I)-I,inv(A-I)*(A+I)+I)
if d==1
    disp('   verifying regularity via verregsing ...')
end
reg=verregsing(B,t);
if reg==1 % outer matrix regular, A is verified to be a P-matrix
    pm=1;
    if d==1
        disp('   verified')
    end
    setround(gr); return % with verified P-property
end
pm=-1; % no verified result
if d==1
    disp('   not verified')
end
setround(gr);
% end of verpmat
%
