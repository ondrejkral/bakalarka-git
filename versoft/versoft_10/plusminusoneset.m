function Y=plusminusoneset(n)
%    PLUSMINUSONESET   Matrix of plus/minus-one vectors.
%
%    Y=plusminusoneset(n) produces a 2^n-by-n matrix Y whose rows are 
%    all plus/minus-one vectors in R^n, with each two adjacent
%    rows of Y differing in exactly one entry.
%
%    COMMENT. This algorithm is employed as a subroutine in full search algorithms
%    that require to perform some operation for all plus/minus-one vectors. 

%    Copyright J. Rohn, 2005.
%
%    Based on the algorithm "ynset" in
%    J. Rohn, A handbook of results on interval linear problems,
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
z=zeros(1,n); 
y=ones(1,n); 
Y=y;
while any(z~=ones(1,n))
    k=find(z==0,1);
    z(1:(k-1))=zeros(1,k-1);
    z(k)=1;
    y(k)=-y(k);
    Y=[Y; y];
end
