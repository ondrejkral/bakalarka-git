function x=verenclinthull(A,b,d) % uses vertex
%    VERENCLINTHULL        Verified enclosure of the interval hull of the solution set
%                          of a system of interval linear equations.
%
%    For a square interval matrix A and a matching interval vector b,
%        x=verenclinthull(A,b,d)
%    either computes a verified enclosure x of the interval hull of the
%    solution set of A*X=b, or yields no verified result. "d" is a search
%    depth parameter; it should be an integer between 0 and 6.
%
%    There are two possible outputs:
%
%    ~isnan(x) :    x is a verified enclosure of the interval hull of A*X=b,
%     isnan(x) :    no verified output.
%
%    Generally said, the command
%        verenclinthull(A,b,d+1)
%    should give a better enclosure (although it is not theoretically proved) than
%        verenclinthull(A,b,d),
%    but it takes approximately twice as much computation time. The initial enclosure
%        verenclinthull(A,b,0)
%    is never worse than that computed by
%        verifylss(A,b).    
%
%    EXAMPLE (Barth and Nuding 1974). Computation done for d=0,1,2,3,4.
%    Observe how the enclosure is successively shrinking down to the optimal one:
% 
%    intval A = 
%    [    2.0000,    4.0000] [   -2.0000,    1.0000] 
%    [   -1.0000,    2.0000] [    2.0000,    4.0000] 
%    intval b = 
%    [   -2.0000,    2.0000] 
%    [   -2.0000,    2.0000] 
%
%    tic, x0=verenclinthull(A,b,0), toc
%    intval x0 = 
%    [  -14.0001,   14.0001] 
%    [  -14.0001,   14.0001] 
%    Elapsed time is 0.137479 seconds.
%
%    tic, x1=verenclinthull(A,b,1), toc
%    intval x1 = 
%    [  -14.0001,   14.0001] 
%    [  -14.0001,   14.0001] 
%    Elapsed time is 0.242569 seconds.
%
%    tic, x2=verenclinthull(A,b,2), toc
%    intval x2 = 
%    [  -12.0953,   12.0953] 
%    [  -12.0953,   12.0953] 
%    Elapsed time is 0.544185 seconds.
%
%    tic, x3=verenclinthull(A,b,3), toc
%    intval x3 = 
%    [   -4.0001,    4.0001] 
%    [   -4.4445,    4.4445] 
%    Elapsed time is 1.113682 seconds.
%
%    tic, x4=verenclinthull(A,b,4), toc
%    intval x4 = 
%    [   -4.0001,    4.0001] 
%    [   -4.0001,    4.0001] 
%    Elapsed time is 2.151929 seconds.
%    COMMENT. Optimal enclosure (interval hull) reached. This, however, is
%    not guaranteed in the general case.
%
%    See also VERINTERVALHULL, VERIFYLSS.

%    Copyright 2005-2008 Jiri Rohn
%
%    The branching mechanism (used in the subfunction VERTEX) is due to S. P. Shary. 
%    The terminal enclosures are computed using the Hansen-Bliek-Rohn
%    enclosure as described in M. Fiedler, J. Nedoma, J. Ramik, J. Rohn 
%    and K. Zimmermann, Linear Optimization Problems with Inexact Data,
%    Springer-Verlag, New York 2006, pp. 72-73.  
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
%    2005-11-03   initial MATLAB version
%    2006-02-15   INTLAB version 
%    2008-12-04   revisited; VERTEX created, embedded
%    2008-12-11   final version of EA, p-coded
%    2008-12-15   final version (html copy)
%
gr=getround;
setround(0);
b=b(:); 
[m,n]=size(A); 
x=repmat(infsup(NaN,NaN),n,1); 
if ~(m==n&&n==length(b)) % wrong data
    return 
end
if ~(nargin==3&&d==round(d)&&0<=d&&d<=6) % d not assigned or not integer or out of the bounds
    d=0; % default
end
if ~isintval(A)
    A=infsup(A,A); % allows for real input
end 
if ~isintval(b)
    b=infsup(b,b); 
end
x=vertex(A,b,0,d); % till depth d
setround(gr);
% end of verenclinthull
%
%%%%%%%%%%%%%%%%%%%%
% Subfunction vertex
%%%%%%%%%%%%%%%%%%%%
%
function x=vertex(A,b,j,d)
% recursive procedure; calls itself
% j current depth (start with j=0)
% d maximal depth allowed
% interval hull guaranteed to be achieved for d=2^(n(n+1))
if j<d
    [m,n]=size(A);
    x=repmat(intval(NaN),n,1);
    x0=ea(A,b);
    Ab=[A b];
    if any(any(Ab.inf<Ab.sup))
        D=Ab.sup-Ab.inf;
        mx=max(max(D)); [i1,j1]=find(D==mx,1);
        if j1<=n
            A1=A; A1(i1,j1)=A.inf(i1,j1);
            x1=vertex(A1,b,j+1,d);
            A2=A; A2(i1,j1)=A.sup(i1,j1);
            x2=vertex(A2,b,j+1,d);
        else
            b1=b; b1(i1)=b.inf(i1);
            x1=vertex(A,b1,j+1,d);
            b2=b; b2(i1)=b.sup(i1);
            x2=vertex(A,b2,j+1,d);
        end
        if ~isnan(x1.inf(1))&&~isnan(x2.inf(1))
            x=hull(x1,x2);
            if ~isnan(x0.inf(1))
                x=intersect(x,x0);
                return
            end
        else
            if ~isnan(x0.inf(1))
                x=x0;
                return
            end
        end
    else
        x=verifylss(A,b);
    end
else % j>=d
    x=ea(A,b);
end
% end of vertex
%