function [x,y]=verabsvaleqn(A,B,b)
%    VERABSVALEQN       Verified solution of the equation A*x+B*abs(x)=b.    
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For real square matrices A, B of the same size and a matching vector b,
%        [x,y]=verabsvaleqn(A,B,b)
%    either produces an interval vector x verified to contain a solution
%    of the equation
%        A*x+B*abs(x)=b,            (1)
%    or finds a real (noninterval) nonzero vector y verified to satisfy
%        abs(A*y)<=abs(B)*abs(y),   (2)
%    or fails (yields vectors x, y of NaN's).
%
%    Possible outcomes:
%
%    ~isnan(x.inf(1))                 :  x is verified to contain a solution of (1),  
%    ~isnan(y(1)))                    :  y is verified to satisfy (2),  
%     isnan(x.inf(1)) && isnan(y(1))) :  no verified result.
%
%    COMMENT. The algorithm is surprisingly effective considering the fact
%    that both (1) and (2) are NP-hard to solve.
%
%    EXAMPLES. Both possible outcomes ale illustrated on two examples with
%    random 100x100 matrices. Due to the length of the resulting vectors,
%    only the first entry is always displayed:
%
%    >> tic, rand('state',4), n=100; A=2*rand(n,n)-1; B=0.1*(2*rand(n,n)-1); b=2*rand(n,1)-1; [x,y]=verabsvaleqn(A,B,b); x(1), y(1), toc
%    intval ans = 
%    [    2.2333,    2.2334] 
%    ans =
%       NaN
%    Elapsed time is 0.246395 seconds.
%    % Solution of (1) found.
%
%    >> tic, rand('state',6), n=100; A=2*rand(n,n)-1; B=0.1*(2*rand(n,n)-1); b=2*rand(n,1)-1; [x,y]=verabsvaleqn(A,B,b); x(1), y(1), toc
%    intval ans = 
%    [       NaN,       NaN] 
%    ans =
%       -0.1947
%    Elapsed time is 0.092539 seconds.
%    % Solution of (2) found.
%
%    See also VERIFYLSS.

%    Copyright 2005-2008 Jiri Rohn.
%
%    Based on an improvement of the algorithm "signaccord" described in
%    J. Rohn, A Handbook of Results on Interval Linear Problems,
%    posted at http://www.cs.cas.cz/~rohn
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
%    2005, Spring  started
%    2005-11-17    absvaleqn (unverified), final version
%    2007-03-26    function absvalverifn
%    2007, Autumn  versions of verabsvaleqn
%    2008, Spring  absvaleqn renamed as ek; outputs changed
%    2008-05-30    version for posting
%
gr=getround;
setround(0);
b=b(:); 
[m,n]=size(A);
% setting default outputs
x=repmat(infsup(NaN,NaN),n,1); 
y=repmat(NaN,n,1); 
if ~(isequal(m,n)&&isequal(size(B),[n,n])&&isequal(length(b),n)&&isreal(A)&&isreal(B)&&isreal(b)&&~isintval(A)&&~isintval(B)&&~isintval(b))
    setround(gr); return % wrong data
end
% finding unverified solution or unverified Oettli-Prager vector
[xunv,yunv]=ek(A,B,b);
% case A: unverified solution found
if ~isnan(xunv(1)) 
    % verification of solution
    xx=absvalverifn(A,B,b,xunv); % verified solution, or an interval vector of NaN's
    if ~isnan(xx.inf(1))
        x=xx;                    % x is a verified solution
        setround(gr); return
    end
end
% case B: unverified Oettli-Prager vector found
if ~isnan(yunv(1))
    % finding a verified Oettli-Prager vector
    yy=infsup(yunv,yunv);
    left =abs(A*yy);             % left-hand side
    right=abs(B)*abs(yy);        % right-hand side
    if all(left.sup<=right.inf)  % Oettli-Prager inequality satisfied
        y=yunv;                  % y is a verified Oettli-Prager vector
        setround(gr); return
    end
end
setround(gr);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions: absvalverifn, ek (outside)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% absvalverifn of 2007-03-26; small changes 2007-09/10, 2008-05-17
function xx=absvalverifn(A,B,b,x)
%    ABSVALVERIFN     Verification of a solution of the equation A*x+B*abs(x)=b.
%
%    For an approximate real solution x of the equation
%        A*x+B*abs(x)=b,     (1)
%    xx=absvalverifn(A,B,b,x)
%    either produces a tight interval vector xx verified to contain a solution 
%    of (1), or fails (yields an interval vector xx of NaN's).

%    Copyright 2007 Jiri Rohn
%
n=size(A,1);
xx=repmat(infsup(NaN,NaN),n,1);
% converting data from real to intval
if ~isintval(A), A=infsup(A,A); end                    
if ~isintval(B), B=infsup(B,B); end
if ~isintval(b), b=infsup(b,b); end
% A, B, b assumed intval quantities already in the 2007-03-26 version
if  isintval(x), return, end
I=eye(n,n);
z=ones(n,1);
for j=1:n
    if x(j)<0, z(j)=-1; end                           % sign vector of x
end
infinf=infsup(repmat(-Inf,n,1),repmat(Inf,n,1));
x1=verifylss(A+B*diag(z),b);                          % enclosure via verifylss
if ~(~isnan(x1.inf(1))&&all(z.*inf(x1)>=0)&&all(z.*sup(x1)>=0))
    x1=infinf;                                        % failure to produce verified output
end                                                   % otherwise x1 is a verified solution
M=inv(I-abs(inv(A)*B));
if isa(M,'double')                                    % bridging the current gap in verifylss
    M=infsup(M,M); 
end
if (~isnan(M.inf(1,1)))&&(all(all(M.inf>=0)))         % M verified nonnegative; guarantees existence and uniqueness of solution
    rad=M*abs(inv(A)*(A*x+B*abs(x)-b));               % enclosure via |x^*-x|<=M*|inv(A)*residual|
    x2=midrad(x,rad.sup);                             % x2 is a verified solution
else
    x2=infinf;                                        % failure to produce verified output
end
x3=intersect(x1,x2);                                  % intersection of the two enclosures
if ~isequal(x3,infinf)                                % at least one enclosure computed
    xx=x3;                                            % verified enclosure of the solution of Ax+B|x|=b
end
% end of absvalverifn
