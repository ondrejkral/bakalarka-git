function [x,As,C]=verintervalhull(A,b,t) 
%    VERINTERVALHULL     Verified interval hull of the solution set of a system of interval linear equations.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%    
%    x=verintervalhull(A,b) computes a verified interval hull x of the
%    solution set of a system of interval linear equations A*x=b, A square, 
%    or yields a vector of NaN's.
%
%    [x,As]=verintervalhull(A,b) gives either a verified interval hull x,
%    or a very tight ("almost thin") interval matrix As which is a part of
%    A and is verified to contain a singular matrix. If no verified result
%    is achieved, then both x and As consist of NaN's.
%
%    [x,As,C]=verintervalhull(A,b) produces additional information in the
%    structure C:
%        C.xun    is an unverified interval hull, issued in case a verified
%                 one could not be computed,
%        C.Asun   is a thin unverified singular matrix, issued in case a verified 
%                 one could not be computed,
%        C.flag   describes verbally the output; it has one of the forms
%                   'verified interval hull computed              '
%                   'interval matrix verified singular            '
%                   'interval hull computed, result not verified  '
%                   'interval matrix singular, result not verified',
%        C.orth   is the number of orthants inspected by the algorithm,
%        C.iter   is the total number of iterations of the core algorithm,
%        C.D      is a (possibly huge) matrix with rows of the form [z xl' xu'], where
%                   z is a plus/minus-one vector, and
%                   infsup(xl,xu) is a verified (possibly not optimal)
%                     enclosure of the intersection of the solution set
%                     with the z-orthant;
%                 if some plus/minus vector z is not present, then it is
%                 verified that the intersection of the solution set with
%                 the z-orthant is empty,
%        C.Xlower is a matrix such that for each i the lower bound x.inf(i) was
%                 attained as the lower bound of the solution of the system
%                   (mid(A)-diag(y)*rad(A)*diag(z))*x=mid(b)+diag(y)*rad(b),
%                 where [y z] is the ith row of C.Xlower,
%        C.Xupper is a matrix such that for each i the upper bound x.sup(i) was
%                 attained as the upper bound of the solution of the system
%                   (mid(A)-diag(y)*rad(A)*diag(z))*x=mid(b)+diag(y)*rad(b),
%                 where [y z] is the ith row of C.Xupper.
%        C.D, C.Xlower, and C.Xupper are empty if x is empty.
%
%    [x,As,C]=verintervalhull(A,b,1) [i.e., with additional input argument "1"] 
%                 is the same as before, but in course of the computation
%                 it also produces screen output of the form 
%                   "Expected remaining time: ... sec."
%                 This is a useful feature, considering the NP-hardness of
%                 the problem; however, it slows down the actual computation.
%
%    See also VERIFYLSS.

%    Copyright 2007-2008 Jiri Rohn
%
%    Based on the algorithm "hull" described in
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
%    2007-11-30   first posted
%    2008-01-07   t=0 for unassigned t added
%
if (nargin==3)&&isequal(t,1)                % t==1: display remaining time 
    t=1;
else
    t=0;
end
[x,As,C]=intervalhull(A,b,t);               % computation done by INTERVALHULL
[m,n]=size(A); 
if isempty(x)                               % output adapted to INTLAB standards
    x=repmat(infsup(NaN,NaN),n,1);
end
if isempty(As)
    As=repmat(infsup(NaN,NaN),m,n);
end
if isempty(C.xun)
    C.xun=repmat(NaN,n,1);
end
if isempty(C.Asun)
    C.Asun=repmat(NaN,m,n);
end
