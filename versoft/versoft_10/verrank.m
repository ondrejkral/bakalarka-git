function [r1,r2]=verrank(A)
%    VERRANK        Verified bounds on the rank of a rectangular real matrix.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a rectangular real matrix A,
%        [r1,r2]=verrank(A)
%    computes integers r1, r2 verified to satisfy
%        r1 <= rank(A) <= r2,
%    where rank(A) is the EXACT rank of A, or fails (yields NaN's). Thus,
%    if r1==r2, then r1 is the verified rank of A. 
%    
%    See also VERSINGVAL, RANK.

%    Copyright 2008 Jiri Rohn.
%
%    Bounds on the rank computed from verified singular values of A using
%    VERSINGVAL.  
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
%    2008-02-05   first version
%    2008-02-14   version for posting
%    2008-04-20   basis part added
%
gr=getround;
setround(0);
r1=NaN; r2=NaN;
if ~(nargin==1&&nargout<=2&&isreal(A)&&~isintval(A)) % wrong data
    setround(gr); return
end
sigma=versingval(A); % main part % interval vector of singular values
if isnan(sigma.inf(1)) % singular values not computed
    setround(gr); return 
end
% singular values computed
r1=length(find(sigma.inf>0)); % number of positive lower bounds
r2=length(find(sigma.sup>0)); % number of positive upper bounds
% next: added 2008-04-20
if isnan(r1)||(~isnan(r1)&&r1<r2) % exact rank not determined
    B=verbasis(A);
    if isnan(B(1,1)) % basis not computed
        setround(gr); return
    end
    % basis computed
    r1=size(B,2);
    r2=r1;
end
setround(gr);
