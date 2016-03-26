function [flag,x,y,C,E]=verlinprogg(A,b,c) % refers to AT
%    VERLINPROGG       Verified linear programming.  
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a real matrix A (full or sparse) and matching real vectors b, c,
%        [flag,x,y,C,E]=verlinprogg(A,b,c)
%    either computes verified optimal solution x, verified dual optimal solution y
%    and verified optimal value C.h of the linear programming problem
%        min c'*x   subject to   A*x=b, x>=0,
%    or verifies (in)feasibility, or verifies unboundedness, or yields no
%    verified result. The respective outcome is always described verbally
%    in the variable flag.
%
%    Possible outputs:
%
%    flag='verified optimum   ' : x is verified to enclose a primal optimal solution, 
%                                 y is verified to enclose a dual optimal solution, 
%                                 C.h is verified to enclose the optimal value,
%                                 C.B is the set of indices of the basic variables,
%                                 C.N is the set of indices of the nonbasic variables,
%                                 C.crit is the final criterial row,
%    flag='verified unbounded ' : x is verified to enclose a primal feasible solution xo, and
%                                 y is verified to enclose a vector yo such that the objective 
%                                   tends to -Inf along the feasible half-line {xo+t*yo | t>=0},
%                                 C fields consist of NaN's,
%    flag='verified feasible  ' : x is verified to enclose a primal feasible solution
%                                   (optimality nor unboundedness could be verified),
%                                 y and C are of NaN's,
%    flag='verified infeasible' : y is verified to enclose a Farkas vector yo satisfying A'*yo>=0, b'*yo<0
%                                   (whose existence proves primal infeasibility),
%                                 x and C are of NaN's,
%    flag='no verified result ' : x, y and C are all of NaN's
%                                   (no verified result could be found),
%    flag='sizes do not match ' : x, y and C are all of NaN's
%                                   (sizes of A, b, c do not match).
%
%    The structured array E explains reasons for NaN output. It has three fields: 
%    E.error, E.where, E.value. 
%
%    See also LINPROG, VERIFYLSS.

%    Copyright 2007-2008 Jiri Rohn
%
%    A part of this work (nondegenerate solutions) was supported by the
%    Czech Republic National Research Program "Information Society",
%    project 1ET400300415. The second part (degenerate solutions) was done
%    during author's employment at the Anglo-American University in Prague,
%    Czech Republic. 
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
%    2007, Spring first working version created
%    2007-05-31   first posted as VERLINPROG; handles nondegenerate solutions only
%    2008-09-10   started again
%    2008-09-11   degenerate solutions handled via new VEROPT using VERSOL
%    2008-10-08   output structured array C created, Xopt added into it (subfunction optsolset)
%    2008-10-11   optsolset included into main file
%    2008-11-08   output of C.Xopt suppressed
%    2008-11-09   C's fields explained in the help
%    2008-11-28   input variable d and all the messages added
%    2008-11-29   AT called by vverlinprog; AT p-coded
%    2008-12-16   d dropped; renamed as verlinprogg; final version (html)  
%
[flag,x,y,C,E]=at(A,b,c); % computation done by AT