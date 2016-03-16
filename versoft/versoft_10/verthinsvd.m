function [U,S,V,E]=verthinsvd(A)
%    VERTHINSVD        Verified thin singular value decomposition of a complex (or real) matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For an m-by-n complex (or real) matrix A, m>=n,
%        [U,S,V,E]=verthinsvd(A)
%    computes (generally complex) m-by-n interval matrix U, a real diagonal
%    n-by-n interval matrix S and an n-by-n interval matrix V that are verified
%    to contain matrices Uo, So, Vo satisfying (in exact arithmetic):  
%        A=Uo*So*Vo',                                                         
%        Uo'*Uo=eye(n,n),
%        Vo'*Vo=eye(n,n),
%        So has nonnegative diagonal entries ordered in nonincreasing order.
%    Hence, Uo, So and Vo form a thin singular value decomposition (SVD) of A. 
%    If A is real, then U and V are real. For s=diag(S), both s.inf and s.sup 
%    are nonnegative and ordered in nonincreasing order. If no verified
%    output is given, then U, S and V consist of NaN's. 
%
%    If m<n, then the decomposition is computed by
%        [U1,S1,V1]=verthinsvd(A');
%        U=V1; S=S1'; V=U1;
%    so that U, S are m-by-m and V is n-by-m and the above properties again
%    hold, this time with
%        Uo'*Uo=eye(m,m),
%        Vo'*Vo=eye(m,m).
%
%    The structure E explains reasons for NaN output. It has three fields:
%    E.error, E.where, E.value.  
%
%    Built-in function.

%    Copyright 2008 Jiri Rohn.
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
%    2008-02-03   first version
%    2008-02-16   version for posting (as versvd)
%    2008-03-28   thin version, verified thin decomposition 
%    2008-03-31   output variable E added; version for posting
%    2008-04-04   transmission of error information from OL via Eol
%    2008-05-30   JK, p-coded, called by VERTHINSVD
%
[U,S,V,E]=jk(A); % computation done by JK