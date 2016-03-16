function iv = ilspencmono( A, b, p, parameter)
%ILSPENCMONO Method for solving interval linear system with dependencies
%based on monotonicity approach.
%   A template for handling this problem is presented. Various attitudes
%   in each step are possible.

% First, we need to compute initial encolsure x*. 
% We can use some of our other methods.
x = ilspenchbr(A,b,p);

% Checking derivative signs of solution with respect to parameters.
D = intval(zeros(length(x),length(p)));
for k = 1:length(p)
    % Setting right side of system to bk - Ak*x, which is constant.
    % Our methods behave differently, depending on data model used.
    % Thus we need to create right side vector with respect to it.
    db = ilspencmakeb(ilspencgetbk(b,k) - ilspencgetak(A,k)*x, 0, p);
    
    % Solving derivate signs system by some of our other methods
    % as proposed in MCM method (Skalna 2008).
    D(:,k) = ilspenchbr(A,db,p);
end

% Now we can comopute tight enclosures for dimensions of x in which 
% is x with respect to p monotonous.
% Several approaches to handling non-monotonous parameters can be applied.
iv = intval(zeros(length(x),1));
switch(parameter)
    case 'NOIMPROVE'
        % no improvement for non-monotonous components
        for i= 1:length(x)
            newx = ilspencmonogetbound(A, b, p, D(i,:));
            iv(i) = newx(i);
        end
    otherwise
        % rise error
end
end


