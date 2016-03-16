function [fcr,E]=verfullcolrank(A)
%    VERFULLCOLRANK      Verified full column rank of a rectangular real matrix. 
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a rectangular real matrix A, 
%        fcr=verfullcolrank(A)
%    either verifies full column rank of A, or verifies column rank
%    deficiency of A, or fails (yields no verified result): 
%
%    fcr= 1      A is verified to have full column rank (i.e., columns of A are linearly independent),
%    fcr= 0      A is verified not to have full column rank (i.e., columns of A are linearly dependent),
%    fcr=-1      no verified result.
%
%    For full row rank, apply VERFULLCOLRANK to A'.
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
%    2007-11-02   first version
%    2008-03-12   version for posting
%    2008-04-05   column rank deficiency part added; E added
%    2008-05-30   ZD, p-coded, called by VERFULLCOLRANK
%
[fcr,E]=zd(A); % computation done by ZD