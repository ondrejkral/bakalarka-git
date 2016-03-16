function [x,E,C]=verhullpatt(A,b) % uses matpatt, jz and, inside it, ea
%    VERHULLPATT     Verified enclosure of the solution set of a system of interval linear equations
%                    with equality of coefficients being taken into account.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
% 
%    For a square interval matrix A and a matching interval vector b,
%        [X,E,C]=verhullpatt(A,b)
%    computes an interval vector X verified to enclose the set
%        { x | Ao*x = bo, Ao parametrized in A, bo in b }    
%    where a matrix Ao in A is called parametrized if it possesses the
%    following property: Ao(i,j)=Ao(k,m) whenever A(i,j)=A(k,m) (that
%    means, the entries of Ao are enforced to be the same whenever their
%    bounds are the same). If no verified enclosure is found, then X consists  
%    of NaN's. The structured array E explains reasons for NaN output. It
%    has three fields: E.error, E.where, E.value. C records some
%    intermediate values. 
%
%    See also VERHULLPARAM, VERINTERVALHULL.

%    Copyright 2008 Jiri Rohn
%
%    The problem is internally reformulated as the problem of solving a
%    system of parametric interval linear equations and then solved as such
%    by VERHULLPARAM.
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
%    2008-11-17   matpatt started
%    2008-11-18   first working version
%    2008-12-31   final version (html)
%
gr=getround;
setround(0);
if nargin==2
    amas=-1; % default: a and -a should be parametrized by different parameters % changed 2008-12-30
  % amas= 1 would mean that a and -a are paramatrized by a single parameter (default amas=-1 set for simplicity)          
end
b=b(:);
[m,n]=size(A);
x=repmat(infsup(NaN,NaN),n,1); % default
E.error='verhullpatt: none';
E.where='NaN';
E.value='NaN';
C.M=repmat(NaN,m,n);
C.Z=NaN; % should be px1, p not determined yet
C.cellarray=NaN; % should be of length 2*(p+1)+1, p not determined yet
C.DxDt=repmat(intval(NaN),1,1); % should be nxp; n, p not determined yet
C.D=repmat(NaN,1,1); % should be nxp; n, p not determined  yet
if ~(m==n&&length(b)==n&&isreal(A)&&isintval(A)&&isreal(b)&&isintval(b))
    E.error='verhullpatt: improper data';
    setround(gr); return
end
% establishing the matrix pattern
[M,Z]=matpatt(A,amas);
C.M=M; % for output
C.Z=Z; % for output
% creating the cell array
p=length(Z);
cellarray=cell(1,2*(p+1)+1); % preallocation
if ~isequal(amas,1) % a and -a regarded as different
    % matrices A0, A1, ..., Ap
    cellarray{1}=zeros(n,n);
    for j=2:p+1
        cellarray{j}=(M==j-1);
    end
    % vectors b0, b1, ..., bp
    cellarray{p+2}=b;
    for j=p+3:2*(p+1);
        cellarray{j}=zeros(n,1);
    end
    % parameter vector t
    cellarray{2*(p+1)+1}=Z;
else % isequal(amas,1) % a and -a parametrized by a single parameter % in the current setting this branch excluded because of amas==-1
    % matrices A0, A1, ..., Ap
    cellarray{1}=zeros(n,n);
    for j=2:p+1
        cellarray{j}=(M==j-1)-(M==-(j-1));
    end
    % vectors b0, b1, ..., bp
    cellarray{p+2}=b;
    for j=p+3:2*(p+1);
        cellarray{j}=zeros(n,1);
    end
    % parameter vector t
    cellarray{2*(p+1)+1}=Z;
end
C.cellarray=cellarray;  % for output
[x,E,CC]=jz(cellarray); % main computation done by JZ
C.xenc=CC.xenc; % for output
C.DxDt=CC.DxDt; % for output
C.D=CC.D;       % for output  
setround(gr);
% end of verhullpatt
%
%%%%%%%%%%%%%%%%%%%%%
% Subfunction matpatt
%%%%%%%%%%%%%%%%%%%%%
%
function [M,Z]=matpatt(A,amas)
% A interval or real
% matrix pattern M, parameters Z
% A(i,j)=Z(M(i,j)) for each i,j
% amas= 1 if a and -a should be parametrized by a single parameter,
% amas~=1 if they should be regarded as different
% 2008-11-17   created
% 2008-11-18   -a added, amas added, real input enabled
ir=0;
if ~isintval(A)
    A=intval(A);
    ir=1; % original A was real
end
if nargin==1
    amas=1; % default: a and -a should be parametrized by a single parameter
end
[m,n]=size(A);
M=zeros(m,n);
Z=intval([]); % Z=[] does not work for matrix A of type intval
p=0;
for i=1:m
    for j=1:n
        if M(i,j)==0  % not yet entered
            a=A(i,j); % current parameter
            p=p+1;    % new length of Z
            Z(p)=a;   % new parameter
            M(i,j)=p; % new pattern record
            for k=i:m
                for l=1:n
                    if isequal(A(k,l),a) % current parameter met
                        M(k,l)=p;        % new pattern record
                    end
                    if isequal(amas,1)   % a and -a parametrized by a single parameter
                        if isequal(A(k,l),-a)&&-a~=a % minus current parameter met
                            M(k,l)=-p;               % new pattern record
                        % else a and -a are regarded as different parameters    
                        end
                    end
                end % of l
            end % of k
        end % of main if
    end % of j
end % of i
Z=Z'; % column vector
if ir==1 % original A was real
    Z=Z.inf;
end
% end of matpatt
%