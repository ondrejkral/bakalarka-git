function [evl,As]=vereigval(A,lambda,t)
%    VEREIGVAL      Verified real eigenvalue of an interval matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a square interval matrix A and a REAL number lambda,
%    [evl,As]=vereigval(A,lambda)
%    verifies whether lambda is an eigenvalue of some matrix in A,
%    or yields no verified result (unfortunately, complex eigenvalues
%    cannot be handled yet): 
%
%    evl= 1           lambda is verified to be an eigenvalue of some matrix in A; 
%                     As is a very tight ("almost thin") interval matrix
%                     which is a part of A and is verified to contain a
%                     real matrix having lambda as an eigenvalue, 
%    evl= 0           lambda is verified not to be an eigenvalue of any matrix in A;
%                     As consists of NaN's,
%    evl=-1           no verified result (data may be wrong).
%
%    [evl,As]=vereigval(A,lambda,1) [i.e., with additional input argument "1"]
%                     is the same as before, but it also produces a screen
%                     output of the form "Expected remaining time: ... sec."
%                     This feature, however, slows down the actual computation.
%
%    Computational experience shows that verifying that lambda is an
%    eigenvalue of some matrix in A (evl=1) is usually very fast; however,
%    verifying that it is not so (evl=0) occasionally may last long since
%    the problem is NP-hard. 
%
%    See also VERREGSING, VEREIGVEC.

%    Copyright 2007 Jiri Rohn
%
%    Based on the section "Real eigenvalues" in
%    J. Rohn, A handbook of results on interval linear problems,
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
gr=getround;
setround(0);
[m,n]=size(A);
evl=-1; As=repmat(infsup(NaN,NaN),m,n);
% checking data
if (m~=n)||~isreal(A)||~isreal(lambda)                % error('wrong data')
    setround(gr); return
end
if ~isintval(A), 
    A=infsup(A,A);                                    % allows for real input 
end
% time display
if (nargin==3)&&isequal(t,1)                          % t==1: display remaining time
    time=1;
else
    time=0;
end
% checking singularity of A-lambda*I
A1=A;
for i=1:n
   A1(i,i)=A1(i,i)-lambda;                            % A1=A-lambda*I    
end
[reg,AAs]=verregsing(A1,time);                        % main part: singularity check of A-lambda*I
% output cases
switch reg
    case  1 
        evl=0;                                        % is not an eigenvalue
    case  0                                           % AAs is verified singular
        for i=1:n
            AAs(i,i)=AAs(i,i)+lambda;                 % AAs+lambda*I
        end
        if in(AAs,A)                                  % AAs part of A
            evl=1; As=AAs;                            % is an eigenvalue; Rohn, SIMAX 1993, Lemma 3.1
        else
            evl=-1;                                   % no verified result
        end
    case -1 
        evl=-1;                                       % no verified result 
end
setround(gr);
