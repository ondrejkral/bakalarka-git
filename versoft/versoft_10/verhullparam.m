function [x,E,C]=verhullparam(varargin) % uses jz and, inside it, ea
%    VERHULLPARAM     Verified enclosure of the solution set of a family of parametric interval linear equations.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
% 
%    For real n-by-n matrices A0, A1, ..., Ap, real n-by-1 vectors b0, b1, ..., bp,
%    and a p-by-1 interval vector t,
%        [X,E,C]=verhullparam(A0,A1,...,Ap,b0,b1,...,bp,t)
%    computes an interval vector X verified to enclose the set
%        { x | A(to)*x = b(to) for some to in t }    (1)
%    where
%        A(to)=A0+to(1)*A1+...+to(p)*Ap
%        b(to)=b0+to(1)*b1+...+to(p)*bp
%    for to in t (hence, (1) is the solution set of a family of parametric
%    interval linear equations). If no verified enclosure is found, then X
%    consists of NaN's. The structured array E explains reasons for NaN output. 
%    It has three fields: E.error, E.where, E.value. The output variable C
%    records some intermediate values. 
%
%    It is also possible to issue the command in the form
%        [X,E,C]=verhullparam(cellarray)
%    where "cellarray" is a cell array satisfying 
%        cellarray{1}=    A0, cellarray{2}=  A1, ..., cellarray{p+1}=  Ap,
%        cellarray{p+2}=  b0, cellarray{p+3}=b1, ..., cellarray{2*p+2}=bp,
%        cellarray{2*p+3)=t.
%    This helps to prepare the data in advance, in particular when the
%    number of parameters p is large (as it may happen with special types
%    of matrices).
%
%    See also VERINTERVALHULL, VERHULLPATT.
%
%    Built-in function.

%    Copyright 2008 Jiri Rohn
%
%    Based on ideas sketched in J. Rohn, A Method for Handling Dependent
%    Data in Interval Linear Systems, Technical Report No. 911, Institute
%    of Computer Science, Academy of Sciences of the Czech Republic, 
%    Prague 2004.
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
%    2008-11-16   started (versolver created; up to construction of D; called verintervalhulldep)
%    2008-11-17   upper and lower loops added, input of a cell array enabled, 
%                 renamed as jz, called by verhullpar; working version
%    2008-11-20   output parameter C added
%    2008-12-04   reworded to "parametric interval linear equations"
%    2008-12-20   jz p-coded, final version (html)
%
[x,E,C]=jz(varargin); % computation done by JZ