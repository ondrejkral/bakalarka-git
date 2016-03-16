function [sol,x]=verbasintnpprob(A,t)
%    VERBASINTNPPROB     Verified solution of the basic NP-complete problem of interval computations.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a square real (noninterval) matrix A,
%    [sol,x]=verbasintnpprob(A)
%    either yields a verified real solution x of the system
%        -e <= A*x <= e,        (1)
%        norm(x,1) >= 1,        (2) 
%    where e is the vector of all ones, or verifies that (1), (2) has no
%    solution, or gives no verified result:
%
%    sol= 1           x is a verified solution of (1), (2),
%    sol= 0           it is verified that (1), (2) has no solution; 
%                     x is a vector of NaN's, 
%    sol=-1           no verified result (A may be not square or not real).
%
%    [sol,x]=verbasintnpprob(A,1) [i.e., with additional input argument "1"]
%                     works in the same way as before, but it also produces screen
%                     output of the form "Expected remaining time: ... sec."
%                     This feature, however, slows down the actual computation.
%
%    COMMENT. In [1], Chapters 2 and 3, it is shown that the problem of
%    solvability of (1), (2) is NP-complete (even for nonnegative
%    positive definite rational matrices), and that NP-completeness or
%    NP-hardness of all the basic interval linear problems can be proved using 
%    this fact. In this way the problem (1), (2) can be considered, at least 
%    theoretically, the basic NP-complete problem of interval computations.
%
%    [1] M. Fiedler, J. Nedoma, J. Ramik, J. Rohn and K. Zimmermann, Linear
%    Optimization Problems with Inexact Data, Springer-Verlag, New York 2006. 

%    Copyright 2008 Jiri Rohn
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
sol=-1; x=repmat(NaN,n,1);
% checking compatibility of data
if ~(m==n&&isreal(A)&&~isintval(A)) % error('matrix not square or not real'),
    setround(gr); return
end
if issparse(A)
    A=full(A); % sparse matrices not implemented yet
end
% time display
if (nargin==2)&&isequal(t,1) % t==1: display remaining time
    time=1;
else
    time=0;
end
tol=1e-10;
e=ones(n,1);
AA=midrad(A,ones(n,n));
% main part: application of verregsing
[reg,As]=verregsing(AA,time); % solution exists if and only if AA is singular; Fiedler et al., proof of Thm. 2.33
% output cases
if reg==-1 % no result
    sol=-1;
    setround(gr); return
end
if reg==1 % AA regular, verified no solution
    sol=0;
    setround(gr); return
end
if reg==0 % AA singular, verified that solution exists
    % finding a solution; solutions are Oettli-Prager vectors of AA taken as null vectors of mid(As)
    X=null(mid(As));
    if isempty(X) % no null vector available
        sol=-1;
        setround(gr); return
    end
    for x1=X % null vectors of mid(As) (candidates for solution) inspected
        x2=x1/(norm(x1,1)-tol); % tol used to ensure norm(x2,1)>=1
        xx=infsup(x2,x2);
        if all(in(A*xx,infsup(-e,e)))&&in(norm(xx,1),infsup(1,Inf)) % xx verified to be a solution
            sol=1;
            x=x2; % noninterval output
            setround(gr); return
        end
    end
    sol=-1;
    setround(gr); return
end

           