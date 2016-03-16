function [L,X,E]=vereig(A)
%    VEREIG        Verified eigenvalues and eigenvectors of a complex (or real) matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a square complex (or real) matrix A,
%        [L,X,E]=vereig(A)
%    computes (generally complex) interval matrices L and X, L diagonal,
%    that are verified to contain matrices Lo, Xo satisfying
%        A*Xo=Xo*Lo
%    in exact arithmetic, where diag(Lo) is the vector of ALL eigenvalues of A
%    and Xo is a matrix of corresponding eigenvectors; L, X are enclosures
%    of these quantities. Multiple eigenvalues are taken into account.
%
%    The vector
%        lam=diag(L)
%    has the following additional property: for each i, j, the intervals
%    lam(i) and lam(j) are either identical, or disjoint. Thus, if all of
%    them are disjoint, then each of them contains exactly one eigenvalue
%    of A.
%
%    If A is real and symmetric, then L, X are real. If A is Hermitian, then L is real. 
%    In these cases both lam.inf and lam.sup are ordered in nondecreasing order.
%
%    The structure E explains reasons for NaN output. It has three fields:
%    E.error, E.where, E.value.  
%
%    EXAMPLE 1 (multiple eigenvalues). The following matrix has a six-tuple
%    eigenvalue 2 corresponding to three Jordan blocks of sizes 1, 2 and 3:
%    A =
%       -60     1    42    -3   -10     4
%       133     0   -92     6    23    -8
%      -186     3   128    -9   -30    12
%       252    -4  -171    14    41   -16
%      -310     5   210   -15   -48    20
%       372    -6  -252    18    60   -22
%    >>  [L,X]=vereig(A); format long, lam=diag(L)
%    intval lam = 
%    [   1.99980487061652 -  0.00019256911560i,   2.00019000884616 +  0.00019256911404i] 
%    [   1.99980487061652 -  0.00019256911560i,   2.00019000884616 +  0.00019256911404i] 
%    [   1.99980487061652 -  0.00019256911560i,   2.00019000884616 +  0.00019256911404i] 
%    [   1.99980487061652 -  0.00019256911560i,   2.00019000884616 +  0.00019256911404i] 
%    [   1.99980487061652 -  0.00019256911560i,   2.00019000884616 +  0.00019256911404i] 
%    [   1.99980487061652 -  0.00019256911560i,   2.00019000884616 +  0.00019256911404i] 
%    Enclosures of equal eigenvalues are equal, as explained above. Low
%    accuracy is caused by multiplicity; the imaginary parts cannot be 
%    filtered out by the program here.
%
%    EXAMPLE 2 (symmetric matrix; Hilbert 4x4). Output for a real symmetric
%    matrix is always real:
%    >> A=hilb(4), [L,X]=vereig(A); format long, lam=diag(L), X
%    A =
%      Columns 1 through 3
%       1.000000000000000   0.500000000000000   0.333333333333333
%       0.500000000000000   0.333333333333333   0.250000000000000
%       0.333333333333333   0.250000000000000   0.200000000000000
%       0.250000000000000   0.200000000000000   0.166666666666667
%      Column 4
%       0.250000000000000
%       0.200000000000000
%       0.166666666666667
%       0.142857142857143
%    intval lam = 
%    [   0.00009670230402,   0.00009670230403] 
%    [   0.00673827360576,   0.00673827360577] 
%    [   0.16914122022144,   0.16914122022146] 
%    [   1.50021428005924,   1.50021428005925] 
%    intval X = 
%      Columns 1 through 2
%    [   0.02919332316478,   0.02919332316479] [   0.17918629053545,   0.17918629053546] 
%    [  -0.32871205576320,  -0.32871205576317] [  -0.74191779062846,  -0.74191779062845] 
%    [   0.79141114583312,   0.79141114583313] [   0.10022813694718,   0.10022813694721] 
%    [  -0.51455274999717,  -0.51455274999714] [   0.63828252819360,   0.63828252819363] 
%      Columns 3 through 4
%    [  -0.58207569949724,  -0.58207569949723] [   0.79260829116376,   0.79260829116377] 
%    [   0.37050218506709,   0.37050218506710] [   0.45192312090159,   0.45192312090160] 
%    [   0.50957863450179,   0.50957863450180] [   0.32241639858182,   0.32241639858183] 
%    [   0.51404827222216,   0.51404827222217] [   0.25216116968824,   0.25216116968825] 
%    Observe the high accuracy of the result.
%    
%    See also VERIFYEIG, EIG.
%
%    Built-in function.

%    Copyright 2008 Jiri Rohn.
%
%    Employs the routine VERIFYEIG by Siegfried M. Rump.
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
%    2008-01-31   first version; handles simple eigenvalues only
%    2008-02-01   second version; multiple eigenvalues allowed (via "orderright")
%    2008-02-06   version for posting
%    2008-02-22   posted as ol (p-coded)
%    2008-03-13   help extended, examples added
%    2008-04-04   output variable E added, nothing else changed
%    2008-05-17   renamed as vereig, ol p-coded
%    2008-05-28   version for posting
%
[L,X,E]=ol(A); % computation done by OL
