function [x,all,E]=verlcpall(M,q,k,t) 
%    VERLCPALL     Verified all solutions of a linear complementarity problem.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For a real square matrix M and a matching real vector q,
%        [X,all,E]=verlcpall(M,q)
%    produces a tight interval matrix X verified to contain as columns 
%    all or some of the solutions of the linear complementarity problem
%        max(x,0)=M*max(-x,0)+q    (1)
%    (in TeX: $x^{+}=Mx^{-}+q$), details being explained below.
%
%    Possible outputs:
%
%    all== 1:   it is verified that X contains as columns all solutions 
%               of (1), all of them verified; thus if X is empty, then 
%               it is verified that (1) has no solution,               
%
%    all==-1:   columns of X are verified solutions of (1), but it is not
%               verified that X contains all the solutions.
%
%    The command
%        [X,all,E]=verlcpall(M,q,k)
%    (i.e., with additional integer input argument k) enforces the computation
%    to be stopped after a prescribed number of solutions k has been found;
%    this can essentially reduce the computation time. In this case it is
%    always all==-1. In particular, use   
%        [X,all,E]=verlcpall(M,q,1)
%    if you are interested in one solution only.
%
%    The command
%        [X,all,E]=verlcpall(M,q,k,1)
%    (i.e., with additional input argument "1") performs the same task, but
%    it also produces screen output of the form 
%        Expected remaining time: ... sec.
%        Number of solutions found: ...
%    This feature, however, slows down the actual computation. Use k=Inf if 
%    you do not wish to bound the number of solutions to be found.
%
%    The structured array E explains details of the output. It has three fields: 
%    E.error, E.where, E.value.
%
%    EXAMPLE 1 (matrix 6-by-6; 14 solutions).
%    >> n=6; rand('state',940); M=2*rnd(n,n)-1, q=2*rnd(n,1)-1, [x,all]=verlcpall(M,q)
%    M =
%       -0.7267    0.3215    0.4428   -0.3401    0.0279   -0.0413
%        0.1583    0.9376   -0.8380    0.5187   -0.2159    0.1201
%        0.8271   -0.2549   -0.8813   -0.3767   -0.4145   -0.3428
%       -0.8641   -0.9621    0.2400    0.5504   -0.3624    0.9646
%       -0.7942    0.9878    0.6991   -0.2685   -0.9371    0.8826
%        0.0365   -0.3602    0.7635    0.8513   -0.8853   -0.4667
%    q =
%        0.4392
%       -0.1077
%        0.9563
%        0.3622
%        0.5596
%        0.4114
%    intval x = 
%      Columns 1 through 5
%    [    0.9386,    0.9387] [   -1.0832,   -1.0831] [   -0.5260,   -0.5259] [    0.4761,    0.4762] [   -0.5097,   -0.5096] 
%    [   -0.6988,   -0.6987] [   -0.4880,   -0.4879] [    0.0625,    0.0626] [   -0.1149,   -0.1148] [    0.0598,    0.0599] 
%    [   -0.7893,   -0.7892] [   -1.4139,   -1.4138] [    1.3281,    1.3282] [    0.9270,    0.9271] [    1.2549,    1.2550] 
%    [   -0.2195,   -0.2194] [   -1.2790,   -1.2789] [   -0.1677,   -0.1676] [    0.2517,    0.2518] [   -0.2111,   -0.2110] 
%    [    1.7427,    1.7428] [    0.8263,    0.8264] [    0.0968,    0.0969] [    0.6730,    0.6731] [   -0.1048,   -0.1047] 
%    [    0.9490,    0.9491] [    2.4433,    2.4434] [    0.5732,    0.5733] [    0.3699,    0.3700] [    0.5168,    0.5169] 
%      Columns 6 through 10
%    [   -0.8018,   -0.8017] [    0.5393,    0.5394] [    0.5397,    0.5398] [   -0.8696,   -0.8695] [   -0.5861,   -0.5860] 
%    [   -0.3242,   -0.3241] [   -0.3410,   -0.3409] [   -0.3415,   -0.3414] [   -0.3802,   -0.3801] [   -0.0316,   -0.0315] 
%    [   -0.9472,   -0.9471] [   -0.2909,   -0.2908] [   -0.2879,   -0.2878] [   -0.6665,   -0.6664] [    1.1487,    1.1488] 
%    [   -1.1696,   -1.1695] [   -0.4911,   -0.4910] [   -0.4862,   -0.4861] [   -0.6765,   -0.6764] [   -0.0659,   -0.0658] 
%    [   -0.6309,   -0.6308] [   -1.0329,   -1.0328] [   -1.0383,   -1.0382] [   -1.2099,   -1.2098] [   -0.3870,   -0.3869] 
%    [    1.4841,    1.4842] [    0.0143,    0.0144] [   -0.0061,   -0.0060] [   -0.6854,   -0.6853] [   -0.2889,   -0.2888] 
%      Columns 11 through 14
%    [    0.6969,    0.6970] [   -1.5772,   -1.5771] [   -0.5519,   -0.5518] [    0.4035,    0.4036] 
%    [   -0.3484,   -0.3483] [   -0.7017,   -0.7016] [    0.0907,    0.0908] [   -0.0022,   -0.0021] 
%    [   -0.4560,   -0.4559] [   -1.3336,   -1.3335] [    1.0958,    1.0959] [    0.6541,    0.6542] 
%    [    1.4469,    1.4470] [    1.1955,    1.1956] [    0.7770,    0.7771] [    1.2087,    1.2088] 
%    [    2.4215,    2.4216] [    3.2669,    3.2670] [    0.9372,    0.9373] [    1.3381,    1.3382] 
%    [   -1.3586,   -1.3585] [   -2.6450,   -2.6449] [   -0.9246,   -0.9245] [   -0.8798,   -0.8797] 
%    all =
%         1
%    Comment: All 14 solutions found, all verified.
%
%    WARNING. The algorithm may generally take up to 2^n steps, where n=size(M,1);
%    it is therefore recommended to run the file for moderate values of n
%    (say, n<=17) only. The exception is the case of a unique solution
%    which, in contrast, can be usually found very quickly. Also, the
%    linear complementarity problem (1) may have up to 2^n solutions, as
%    it can be shown by the simple example 
%        max(x,0)=-eye(n,n)*max(-x,0)+ones(n,1),
%    whose solutions are all the plus/minus vectors in R^n, i.e., 2^n of
%    them.
%
%    EXAMPLE 2 (matrix 6-by-6, 64 solutions).
%    >> n=6; rand('state',3); M=-eye(n,n)+0.1*(2*rand(n,n)-1), q=rand(n,1), [x,all]=verlcpall(M,q); sols=size(x,2), all
%    M =
%       -0.9968    0.0643   -0.0134    0.0826    0.0476    0.0795
%       -0.0550   -1.0261    0.0222    0.0181    0.0023   -0.0633
%       -0.0633   -0.0941   -1.0903    0.0150   -0.0610   -0.0957
%       -0.0567   -0.0616    0.0615   -1.0999    0.0735   -0.0818
%       -0.0146   -0.0506    0.0017   -0.0682   -0.9732    0.0043
%        0.0941    0.0134   -0.0369   -0.0450   -0.0210   -1.0332
%    q =
%        0.4341
%        0.2847
%        0.6537
%        0.1906
%        0.0641
%        0.2679
%    sols =
%        64
%    all =
%         1
%    Comment: The equation has exactly 64 (=2^6) solutions, all of them found verified.
%
%    See also VERABSVALEQNALL. 

%    Copyright Jiri Rohn, 2008.
%
%    Solves an equivalent absolute value equation (I+M)*x+(I-M)*abs(x)=2*q 
%    by VERABSVALEQNALL.
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
%    2008-11-02   first version
%    2008-11-05   input variable k and output variable E added, help written
%    2008-11-23   examples added
%    2008-11-24   final version (html)
%
gr=getround;
setround(0);
% setting default outputs
x=[];
all=-1; 
E.error='verlcpall: none';
E.where='NaN';
E.value='NaN';
if (nargin<4)||((nargin==4)&&~isequal(t,1))
    t=0; % no time display
end
if nargin<3
    k=Inf; % default: all solutions
end
if nargin<2
    E.error='verlcpall: few input data';
    setround(gr); return
end
q=q(:); 
[m,n]=size(M);
if ~(isequal(m,n)&&isequal(length(q),n))
    E.error='verlcpall: sizes do not match';
    setround(gr); return
end
I=eye(n,n); I=infsup(I,I); 
q=infsup(q,q);
% rewriting x^{+}=M*x^{-}+q as (I+M)*x+(I-M)*|x|=2*q, to be solved by verabsvaleqnall
% interval data for verabsvaleqnall
A=I+M; 
B=I-M;
b=2*q;
[x,all,Everabsvaleqnall]=verabsvaleqnall(A,B,b,k,t);
E=Everabsvaleqnall;
E.error=['verlcpall:' E.error];
% end of verlcpall
%