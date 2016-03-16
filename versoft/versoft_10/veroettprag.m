function [oe,Ao,bo]=veroettprag(A,b,x)
%    VEROETTPRAG    Verified solution of the Oettli-Prager inequality, and the data to it.
% 
%    For a rectangular interval matrix A and matching interval vector b and
%    real vector x,
%    [oe,Ao,bo]=veroettprag(A,b,x)
%    either verifies that x solves the Oettli-Prager inequality
%        abs(A.mid*x-b.mid) <= A.rad*abs(x)+b.rad,      (1)                                            
%    or verifies that it is not so, or yields no verified result.
%    In the positive case it also finds data Ao, bo to x. 
%
%    Possible outputs:
%
%    oe= 1       x is verified to satisfy (1); 
%                Ao is a very tight interval matrix which is a part of A and
%                bo is a very tight interval vector which is a part of b
%                such that x is verified to satisfy A1*x=b1 (in exact arithmetic) 
%                for some A1 in Ao, b1 in bo,
%    oe= 0       x is verified not to satisfy (1); Ao, bo consist of NaN's,
%    oe=-1       no verified result.
%
%    COMMENT. The Oettli-Prager theorem (Num. Math. 1964) says that a
%    vector x satisfies the inequality (1) if and only if A1*x=b1 holds for
%    some A1 in A, b1 in b. Hence, the main task here is not only to verify
%    (1), but also to find A1, b1. They are enclosed in the computed Ao, bo.
%
%    See also VERINTERVALHULL.

%    Copyright 2008 Jiri Rohn
%
%    Based on Theorem 2.9 and Proposition 2.10 in
%    M. Fiedler, J. Nedoma, J. Ramik, J. Rohn and K. Zimmermann, Linear
%    Optimization Problems with Inexact Data, Springer-Verlag, New York
%    2006.
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
b=b(:); 
oe=-1;
Ao=repmat(infsup(NaN,NaN),m,n);
bo=repmat(infsup(NaN,NaN),m,1);
if nargin<3||m~=length(b)||n~=length(x)||~isreal(A)||~isreal(b) 
    setround(gr); return
end
if ~isintval(A) % allows for real input
    A=infsup(A,A); 
end
if ~isintval(b) 
    b=infsup(b,b); 
end
z=ones(n,1); % preallocation
for i=1:n
    if x(i)>=0, z(i)=1; else z(i)=-1; end % sign vector of x
end
xi=infsup(x,x);
[Ac,Delta]=vermidrad(A);
[bc,delta]=vermidrad(b);
oeprl=abs(Ac*xi-bc);               % Oettli-Prager inequality, left  side
oeprr=Delta*abs(xi)+delta;         % Oettli-Prager inequality, right side
if all(oeprl.sup<=oeprr.inf)       % Oettli-Prager inequality satisfied
    y=(Ac*xi-bc)./oeprr;           % componentvise division
    y(find(isnan(y)))=infsup(1,1); % case od both numerator and denominator being zero
    Ao=Ac-(diag(y)*Delta)*diag(z); % construction of Ao
    bo=bc+diag(y)*delta;           % construction of bo 
    Ao=intersect(Ao,A);            % Ao made part of A
    bo=intersect(bo,b);            % bo made part of b 
    if ~any(any(isnan(Ao)))&&~any(any(isnan(bo))) % intersections nowhere empty
        oe=1;                      % x verified to satisfy the inequality, Ao, bo found
        setround(gr); return                                        
    else
        oe=-1;                     % x verified to satisfy the inequality, Ao, bo not found
        Ao=repmat(infsup(NaN,NaN),m,n);
        bo=repmat(infsup(NaN,NaN),m,1);
        setround(gr); return 
    end
end
if any(oeprr.sup<oeprl.inf)       
    oe=0;                          % x verified not to satisfy the inequality
    setround(gr); return
end
setround(gr);                      % no verified result
%
% Subfunction
%
function [Ac,Delta]=vermidrad(A)
% computes verified midpoint and radius of A
% Ac, Delta are intval quantities
if ~isintval(A)
    Ac=infsup(A,A);
    Delta=infsup(zeros(size(A)),zeros(size(A)));
else
    Al=infsup(A.inf,A.inf);
    Au=infsup(A.sup,A.sup);
    Ac=   (Al+Au)/2;        
    Delta=(Au-Al)/2;        
end
