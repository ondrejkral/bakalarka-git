function [x,all,E]=verabsvaleqnall(A,B,b,k,t) % uses performwithy, xysolution, absvalverifn, sgn, timeprint
%    VERABSVALEQNALL     Verified all solutions of the equation A*x+B*abs(x)=b.
%
%    This is an INTLAB file. It requires to have INTLAB installed under
%    MATLAB to function properly.
%
%    For real square matrices A, B of the same size and a matching real vector b,
%        [X,all,E]=verabsvaleqnall(A,B,b)
%    produces a tight interval matrix X verified to contain as columns all
%    or some of the verified solutions of the equation 
%        A*x+B*abs(x)=b   (1)
%    (explained below).
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
%        [X,all,E]=verabsvaleqnall(A,B,b,k)
%    (i.e., with additional integer input argument k) enforces the computation
%    to be stopped after a prescribed number of solutions k has been found;
%    this can essentially reduce the computation time. In this case it is
%    always all==-1. In particular, use   
%        [X,all,E]=verabsvaleqnall(A,B,b,1)
%    if you are interested in one solution only.
%
%    The command
%        [X,all,E]=verabsvaleqnall(A,B,b,k,1) 
%    (i.e., with additional input argument "1") performs the same task, but
%    it also produces screen output of the form 
%        Expected remaining time: ... sec.
%        Number of solutions found: ...
%    This feature, however, slows down the actual computation. Use k=Inf if 
%    you do not wish to bound the number of solutions to be found.
%
%    The structured array E explains details of the output. It has three
%    fields: E.error, E.where, E.value.
%
%    EXAMPLE 1 (matrices 7-by-7: 10 solutions).
%    >> tic, n=7; rand('state',671); A=2*rand(n,n)-1, B=2*rand(n,n)-1, b=2*rand(n,1)-1, [x,all]=verabsvaleqnall(A,B,b), toc
%    A =
%       -0.1479   -0.5985   -0.2265   -0.2292   -0.2426   -0.4978    0.4772
%        0.3503    0.7914   -0.8554    0.2560   -0.4149   -0.3221   -0.5674
%       -0.8144    0.8176   -0.9111   -0.9181    0.1953   -0.9376    0.0201
%        0.1143   -0.8706   -0.1203    0.5198   -0.6242   -0.7633   -0.1536
%        0.7850   -0.7964    0.6195   -0.5218    0.9041    0.7736    0.9708
%       -0.4198   -0.5983    0.9180   -0.5057   -0.6677    0.1967    0.0734
%       -0.1962    0.6255   -0.3860    0.1035    0.4396   -0.7893   -0.9860
%    B =
%       -0.8464   -0.5703   -0.9208   -0.0867    0.2831    0.9318    0.8203
%       -0.7984    0.3861   -0.1074   -0.1288    0.8478    0.8475    0.8466
%        0.3445    0.4156    0.7606   -0.4585    0.9195    0.0428    0.0485
%       -0.1394    0.8962   -0.2990   -0.2622   -0.6214   -0.5709   -0.1978
%        0.8221    0.1798   -0.2713    0.9308   -0.9663    0.9149   -0.0731
%        0.8508   -0.2720   -0.7906   -0.8783    0.5006   -0.9402    0.6437
%        0.7253    0.0865    0.5792   -0.1374   -0.0348    0.4932   -0.2036
%    b =
%       -0.6525
%        0.3719
%        0.6019
%       -0.3199
%        0.2327
%       -0.3168
%        0.5135
%    intval x = 
%      Columns 1 through 5
%    [    0.1483,    0.1484] [   -0.6616,   -0.6615] [   -4.3204,   -4.3203] [   -1.8898,   -1.8897] [    0.2118,    0.2119] 
%    [    0.4041,    0.4042] [    0.5318,    0.5319] [   -0.2405,   -0.2404] [    0.4361,    0.4362] [    0.3699,    0.3700] 
%    [    0.7863,    0.7864] [    0.5816,    0.5817] [   -1.2115,   -1.2114] [   -0.3516,   -0.3515] [    0.1696,    0.1697] 
%    [    0.1354,    0.1355] [    0.9473,    0.9474] [    6.0073,    6.0074] [    2.5083,    2.5084] [    0.0233,    0.0234] 
%    [   -0.1988,   -0.1987] [   -0.3511,   -0.3510] [   -2.2161,   -2.2160] [   -0.7775,   -0.7774] [    0.2023,    0.2024] 
%    [   -0.3220,   -0.3219] [   -0.2792,   -0.2791] [   -0.1361,   -0.1360] [   -0.0891,   -0.0890] [   -0.0745,   -0.0744] 
%    [    0.2679,    0.2680] [    0.6275,    0.6276] [    2.8806,    2.8807] [    1.2929,    1.2930] [    0.0600,    0.0601] 
%   Columns 6 through 10
%    [    0.2842,    0.2843] [   -1.9018,   -1.9017] [    0.2798,    0.2799] [    0.1584,    0.1585] [    0.1048,    0.1049] 
%    [    0.2851,    0.2852] [   -0.3675,   -0.3674] [    0.2884,    0.2885] [    0.3711,    0.3712] [    0.3815,    0.3816] 
%    [   -0.0842,   -0.0841] [   -0.7374,   -0.7373] [   -0.0792,   -0.0791] [    0.1642,    0.1643] [    0.2708,    0.2709] 
%    [   -0.0107,   -0.0106] [    2.2569,    2.2570] [   -0.0202,   -0.0201] [   -0.1676,   -0.1675] [   -0.2814,   -0.2813] 
%    [    0.2235,    0.2236] [   -1.0900,   -1.0899] [    0.2207,    0.2208] [    0.1049,    0.1050] [   -0.0258,   -0.0257] 
%    [    0.0125,    0.0126] [    0.4787,    0.4788] [    0.0113,    0.0114] [   -0.0478,   -0.0477] [   -0.0704,   -0.0703] 
%    [    0.0044,    0.0045] [    0.8552,    0.8553] [   -0.0032,   -0.0031] [   -0.0900,   -0.0899] [   -0.1584,   -0.1583] 
%    all =
%         1
%    Elapsed time is 1.017689 seconds.
%    Comment: All 10 verified solutions have been found.
%
%    WARNING. The algorithm may generally take up to 2^n steps, where n=size(A,1);
%    it is therefore recommended to run the file for moderate values of n
%    (say, n<=17) only. The exception is the case of a unique solution
%    which, in contrast, can be usually found very quickly. Also, the
%    equation (1) may have up to 2^n solutions, as shown by the simple example
%        zeros(n,n)*x+eye(n,n)*abs(x)=ones(n,1),
%    whose solutions are all the plus/minus vectors in R^n, i.e., 2^n of
%    them.
%
%    EXAMPLE 2 (matrices 10-by-10: 1024 solutions).
%    >> tic, n=10; rand('state',1); A=0.1*(2*rand(n,n)-1), B=rand(n,n); B=inv(B), b=rand(n,1), [x,all]=verabsvaleqnall(A,B,b); sols=size(x,2), all, toc
%    A =
%      Columns 1 through 7
%        0.0906    0.0797    0.0537   -0.0211   -0.0119    0.0939    0.0799
%        0.0408   -0.0142   -0.0881    0.0383    0.0994   -0.0393    0.0878
%        0.0908   -0.0601    0.0254    0.0139    0.0125   -0.0252   -0.0452
%        0.0196   -0.0394   -0.0470   -0.0314   -0.0068   -0.0303   -0.0046
%        0.0681    0.0077   -0.0375    0.0190   -0.0162    0.0402   -0.0493
%       -0.0114    0.0820    0.0045   -0.0452   -0.0577    0.0229   -0.0276
%        0.0674    0.0051   -0.0183   -0.0904   -0.0528    0.0518    0.0054
%        0.0037   -0.0386    0.0786    0.0676   -0.0372    0.0598   -0.0238
%       -0.0956   -0.0931    0.0148   -0.0795    0.0740   -0.0569    0.0647
%       -0.0248    0.0431    0.0136    0.0457   -0.0496    0.0787   -0.0371
%      Columns 8 through 10
%        0.0326   -0.0732    0.0081
%        0.0295    0.0870    0.0821
%       -0.0059   -0.0965    0.0162
%        0.0436   -0.0878   -0.0590
%       -0.0962    0.0015    0.0768
%       -0.0043   -0.0246   -0.0429
%        0.0979   -0.0795   -0.0703
%        0.0955    0.0010    0.0952
%       -0.0004   -0.0540   -0.0464
%       -0.0918   -0.0612   -0.0713
%    B =
%      Columns 1 through 7
%       -0.2352   -0.5766    1.8918   -0.4172   -0.0450   -1.4132    0.1931
%       -0.2286    0.1210   -0.1436   -0.6798   -0.0777    0.1745   -0.3110
%       -0.7952   -0.0076    0.8798   -0.4007    0.3096   -0.1172    0.8636
%        0.3066   -0.1892   -0.0645    0.5488    0.9801   -0.5211   -0.3190
%        0.6403    1.4887    0.3371   -0.2049   -0.0406    0.0974    0.1317
%       -0.0278    0.8313   -0.5236    0.7833   -0.5226    0.1913   -0.0099
%        0.5736    0.2986    0.1326   -0.5930   -0.3257   -1.0060    0.6667
%       -0.0876   -1.5485   -1.2575    1.3963    0.2402    1.0701   -1.1621
%       -0.1125   -0.5930   -0.4584    0.3201   -0.3340    0.6841   -0.4610
%        0.3664   -0.4505   -1.2536    0.1229    0.2388    1.2614   -0.0462
%      Columns 8 through 10
%       -0.4376   -0.1736    0.9970
%        0.2743    0.3365    0.5909
%        0.1284   -0.8769    0.1664
%       -0.5032    0.7230   -0.4455
%       -0.6928   -0.7540   -0.5845
%       -0.2025   -0.0502    0.0523
%       -0.1276    0.3920    0.0121
%        1.5229    0.5034   -0.5145
%        0.4155    1.3115   -0.6138
%        0.3281   -0.4432    0.1292
%    b =
%        0.1951
%        0.9186
%        0.7427
%        0.0069
%        0.8383
%        0.3371
%        0.0645
%        0.7455
%        0.3865
%        0.2190
%    sols =
%            1024
%    all =
%         1
%    Elapsed time is 8.318330 seconds.
%    Comment: The equation has exactly 1024 (=2^(10)) solutions, all of them found verified.
%
%    See also VERABSVALEQN, VERLCPALL, VERIFYLSS. 

%    Copyright Jiri Rohn, 2007-2008.
%
%    Looks for verified solutions of A*x+B*abs(x)=b through full search
%    over all systems (A+B*diag(y))*x=b, diag(y)*x>=0, y a plus/minus-one
%    vector (ynset algorithm).
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
%    2007-10 to 2007-12-10  original version (unpublished)
%    2008-11-02   changed to admit interval data (because of verlcp)
%    2008-11-04   output variable E added, case of unique solution handled
%                 via verabsvaleqn, input variable k added, help written 
%    2008-11-05   special call of verabsvaleqn in case of k==1 added
%    2008-11-07   the same elaborated for the case of a call by verlcpall;
%                 print of the number of solutions currently found 
%    2008-11-23   examples added; final version (html)
%
gr=getround;
setround(0);
% setting default outputs
x=[];   % x empty; later to be incremented
all=-1; % all==1: verified that all solutions found
E.error='verabsvaleqnall: none';
E.where='NaN';
E.value='NaN';
if (nargin<5)||((nargin==5)&&~isequal(t,1))
    t=0; % no time display
end
if nargin<4
    k=Inf; % default: all solutions
end
if nargin<3
    E.error='verabsvaleqnall: few input data';
    setround(gr); return
end
b=b(:); 
[m,n]=size(A);
kk=k; % in order input k not to interfere with k in the ynset loop
if ~(isequal(m,n)&&isequal(size(B),[n,n])&&isequal(length(b),n))
    E.error='verabsvaleqnall: sizes do not match';
    setround(gr); return 
end
if nargin>=4&&~((kk>0&&round(kk)==kk)||kk==Inf) % kk nonpositive or not integer, not infinity
    E.error='verabsvaleqnall: wrong bound on the number of solutions';
    setround(gr); return 
end
I=eye(n,n);
% next lines changed against verabsvaleqnall of 2007-12-10 (to admit interval data A, B, b) 
% A1=infsup(A,A); B1=infsup(B,B); % old version
if ~isintval(A)
    A1=infsup(A,A);
else
    A1=A;
end
if ~isintval(B)
    B1=infsup(B,B);
else
    B1=B;     
end
N=inv(A1);
M=inv(I-abs(N*B1));
if isa(M,'double') % bridging the current gap in verifylss
    M=infsup(M,M); 
end
if ~isnan(M.inf(1,1))&&~any(any(M.inf<0))
    reg=1; % verified regular (unique solution)
else
    reg=-1;
end
tic, sh=1; told=-sh; % put separately because of performwithy under t~=1, and of premature stop
if isequal(t,1) % display of remaining time required
    tic, sh=1;  % time shift (used in "timeprint")
    told=-sh;
end
if reg==1||kk==1 % one solution to be outputted
    if ~isintval(A)&&~isintval(B)&&~isintval(b) % all data real
        xr=verabsvaleqn(A,B,b); % verabsvaleqn requires noninterval data
        if ~isnan(xr.inf(1)) % verabsvaleqn successful
            x=xr;
            if reg==1
                all=1; % solution unique
                E.error='verabsvaleqnall: none; full search stopped prematurely, unique solution found';
            else % kk==1
                all=-1;
                E.error='verabsvaleqnall: none; full search stopped prematurely, prescribed number of solutions reached';
            end
            if isequal(t,1)
                clc
                % disp(['Expected remaining time: ', num2str(0), ' sec.'])
                disp(['Elapsed time: ', num2str(round(toc)), ' sec.'])
                disp(['Number of solutions found: ', num2str(size(x,2))])
            end
            setround(gr); return % and exit verabsvaleqnall
        end
    else % not all data real (case of verlcpall)
        % making real input A2, B2, b2
        if ~isintval(A), A2=A; else A2=mid(A); end
        if ~isintval(B), B2=B; else B2=mid(B); end
        if ~isintval(b), b2=b; else b2=mid(b); end
        xr0=verabsvaleqn(A2,B2,b2); % verabsvaleqn requires noninterval data
        if ~isnan(xr0.inf(1)) % solution with noninterval data found
            xr=absvalverifn(A,B,b,xr0,M,N); % verification of the solution xr0 based on noninterval data
            if ~isnan(xr.inf(1)) % verification successful
                x=xr;
                if reg==1
                    all=1; % solution unique
                    E.error='verabsvaleqnall: none; full search stopped prematurely, unique solution found';
                else % kk==1
                    all=-1;
                    E.error='verabsvaleqnall: none; full search stopped prematurely, prescribed number of solutions reached';
                end
                if isequal(t,1)
                    clc
                    % disp(['Expected remaining time: ', num2str(0), ' sec.'])
                    disp(['Elapsed time: ', num2str(round(toc)), ' sec.'])
                    disp(['Number of solutions found: ', num2str(size(x,2))])
                end
                setround(gr); return % and exit verabsvaleqn
            end
        end
    end
end
iter=0; % iteration counter
allfound=1; % allfound will remain 1 if it is verified that all solutions have been found
% ynset algorithm starts
% initial y for ynset
z=zeros(1,n);
x0=verifylss(A,b);
if ~isnan(x0.inf(1))
    x0=mid(x0); % made real to get sign
    y= sgn(x0); % initial y taken as the sign vector of the solution of Ax=b (supported by the result in Fiedler et al.)
else % default setting if verifylss fails
    y=ones(1,n);
end
% perform the task with initial y
[x,iter,all,allfound,told,exit,Eperformwithy]=performwithy(A,B,b,y,M,N,x,reg,t,all,allfound,n,iter,sh,told);
if exit==1
    all=1; % to be sure
    E.error='verabsvaleqnall: none; full search stopped prematurely, unique solution found';
    clc
    disp(['Elapsed time: ', num2str(round(toc)), ' sec.'])
    disp(['Number of solutions found: ', num2str(size(x,2))])
    setround(gr); return
end
if size(x,2)==kk % prescribed bound reached
    all=-1;
    E.error='verabsvaleqnall: none; full search stopped prematurely, prescribed number of solutions reached';
    clc
    disp(['Elapsed time: ', num2str(round(toc)), ' sec.'])
    disp(['Number of solutions found: ', num2str(size(x,2))])
    setround(gr); return
end
if ~isequal(Eperformwithy.where,'NaN')
    E=Eperformwithy;
    E.error=['verabsvaleqnall:' E.error];
end
% ynset algorithm 
while any(z~=ones(1,n))
    % determining new y
    k=find(z==0,1);
    z(1:(k-1))=zeros(1,k-1); 
    z(k)=1;
    y(k)=-y(k);               
    % perform the task with y
    [x,iter,all,allfound,told,exit,Eperformwithy]=performwithy(A,B,b,y,M,N,x,reg,t,all,allfound,n,iter,sh,told);
    if exit==1
        all=1; % to be sure
        E.error='verabsvaleqnall: none; full search stopped prematurely, unique solution found';
        clc
        disp(['Elapsed time: ', num2str(round(toc)), ' sec.'])
        disp(['Number of solutions found: ', num2str(size(x,2))])
        setround(gr); return
    end
    if size(x,2)==kk % prescribed bound reached
        all=-1;
        E.error='verabsvaleqnall: none; full search stopped prematurely, prescribed number of solutions reached';
        clc
        disp(['Elapsed time: ', num2str(round(toc)), ' sec.'])
        disp(['Number of solutions found: ', num2str(size(x,2))])
        setround(gr); return
    end
    if ~isequal(Eperformwithy.where,'NaN')
        E=Eperformwithy;
        E.error=['verabsvaleqnall:' E.error];
    end
end % of ynset
if isequal(allfound,1)
    all=1; % verified that all solutions found
end
if isequal(t,1) % remaining time==0
    clc
    % disp(['Expected remaining time: ', num2str(0), ' sec.'])
    disp(['Elapsed time: ', num2str(round(toc)), ' sec.'])
    disp(['Number of solutions found: ', num2str(size(x,2))])
end
setround(gr);
% full search performed
% end of verabsvaleqnall
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions: performwithy, xysolution, absvalverifn, sgn, timeprint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [x,iter,all,allfound,told,exit,E]=performwithy(A,B,b,y,M,N,x,reg,t,all,allfound,n,iter,sh,told) % uses xysolution
% performs the main task of the ynset algorithm (called twice)
% N=inv(A);
% M=inv(I-abs(N*B));
% x is the matrix whose columns are solutions already found 
% only exit and E are new variables
E.error='performwithy: none';
E.where='NaN';
E.value='NaN';
iter=iter+1; 
[xy,eexist,Exysolution]=xysolution(A,B,b,y,M,N);
if ~isnan(xy.inf(1))
    if isempty(x) 
        x=xy;
        if reg==1 % spectral condition satisfied, unique solution found
            all=1;
            if isequal(t,1)
                % clc
                % disp(['Expected remaining time: ', num2str(0), ' sec.'])
                disp(['Elapsed time: ', num2str(round(toc)), ' sec.'])
                disp(['Number of solutions found: ', num2str(size(x,2))])
            end
            exit=1; return % and exit verabsvaleqnall
        end
    else
        if ~ismember(xy',x','rows') % add only if not there yet
            x=[x xy]; 
        end
    end
else
    E=Exysolution;
    E.error=['performwithy:' E.error];
end
if ~(isequal(eexist,1)||isequal(eexist,0)) % verified (non)existence of xy-solution not established
    allfound=allfound-1;
end
if isequal(t,1) % "Expected remaining time" printed only here
    told=timeprint(2^n-iter,iter,sh,told,x);
end
exit=-1;
% end of performwithy
%
function [xy,exist,E]=xysolution(A,B,b,y,M,N) % uses absvalverifn
% solution of (A+B*diag(y))*x=b, diag(y)*x>=0 for particular y
% M, N passed to absvalverifn
% N=inv(A);
% M=inv(I-abs(N*B));
% xy interval vector (verified or of NaN's) 
% exist=1 for verified solution, =0 for verified nonexistence, =-1 otherwise 
n=length(b);
xy=repmat(infsup(NaN,NaN),n,1); exist=-1; 
E.error='xysolution: none';
E.where='NaN';
E.value='NaN';
xcurr=verifylss(A+B*diag(y),b); % current solution, interval quantity
if isnan(xcurr.inf(1))&&~isintval(A)&&~isintval(B)&&~isintval(b) % verifylss fails
    xcurr=pinv(A+B*diag(y))*b; % replaced by lsq solution; real quantity
    xcurr=absvalverifn(A,B,b,xcurr,M,N); % verification; interval quantity
end
if ~isnan(xcurr.inf(1)) 
    yx=diag(y)*xcurr;
    if all(yx.inf>=0) 
        xy=xcurr; exist=1; % verified solution
        return
    end
    if all(yx.sup>=0) % possibility of solution
        xx=absvalverifn(A,B,b,mid(xcurr),M,N); % verification 
        if ~isnan(xx.inf(1))
            xy=xx; exist=1; % verified solution
            return
        end
    end
    if any(yx.sup<0) % solution verified not to exist
        exist=0; % verified no solution for y
        return
    end
end
if isnan(xcurr.inf(1)) 
    E.error='xysolution: verifylss failed';
    E.where=y;
    E.value=xcurr;
else % ~isnan(xcurr.inf(1)) 
    E.error='xysolution: verification failed';
    E.where=y;
    E.value=xcurr;
end
% end of xysolution
%
function xx=absvalverifn(A,B,b,x,M,N) % uses only verifylss
%    ABSVALVERIFN     Verification of a solution of the equation A*x+B*abs(x)=b.
%
%    For an approximate real solution x of the equation
%        A*x+B*abs(x)=b,     (1)
%    xx=absvalverifn(A,B,b,x)
%    either produces a tight interval vector xx verified to contain a solution of (1), 
%    or fails (yields an interval vector xx of NaN's).
%    M, N passed in order not to be computed anew.
%    N=inv(A);
%    M=inv(I-abs(N*B));

%    Copyright 2007 Jiri Rohn
%
n=size(A,1);
xx=repmat(infsup(NaN,NaN),n,1);
if ~isintval(A), A=infsup(A,A); end                   % converting data from real to intval 
if ~isintval(B), B=infsup(B,B); end
if ~isintval(b), b=infsup(b,b); end
if  isintval(x), return, end
z=sgn(x);
infinf=infsup(repmat(-Inf,n,1),repmat(Inf,n,1));
x1=verifylss(A+B*diag(z),b);                          % enclosure via verifylss
if ~(~isnan(x1.inf(1))&&all(z.*inf(x1)>=0)&&all(z.*sup(x1)>=0))
    x1=infinf;                                        % failure to produce verified output
end
if (~isnan(M.inf(1,1)))&&(all(all(M.inf>=0)))         % M verified nonnegative
    rad=M*abs(N*(A*x+B*abs(x)-b));                    % N=inv(A)
    x2=midrad(x,rad.sup);                             % enclosure via |x^*-x|<=M*|inv(A)*residual|
else
    x2=infinf;                                        % failure to produce verified output
end
x3=intersect(x1,x2);                                  % intersection of the two enclosures
if x3~=infinf                                         % at least one enclosure computed
    xx=x3;                                            % verified enclosure of the solution of Ax+B|x|=b
end
% end of absvalverifn
%
function z=sgn(x)
% signum of x for real or intval x
% x column or row, z column
n=length(x);
z=zeros(n,1);
if ~isintval(x) % real vector
    for i=1:n
        if x(i)>=0
            z(i)=1;
        else
            z(i)=-1;
        end
    end
else % intval vector
    for i=1:n
        if x.inf(i)>0
            z(i)=1;
        end
        if x.sup(i)<0
            z(i)=-1;
        end
        % otherwise z(i)=0
    end
end
% end of sgn
%
function told=timeprint(z,d,sh,told,x)
% prints remaining time
if toc-told>=sh
    clc
    disp(['Expected remaining time: ', num2str(round((toc/d)*z)), ' sec.'])
    disp(['Number of solutions found: ', num2str(size(x,2))])
    told=toc;
end
% end of timeprint
%