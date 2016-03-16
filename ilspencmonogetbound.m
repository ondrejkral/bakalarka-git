function iv = ilspencmonogetbound( A, b, p, derivates)
%ILSPENCMONOGETBOUND Compouting bounds of i-th component of x with respect
%to parameters monotonicity.

lowerp = intval(zeros(length(derivates),1));
upperp = intval(zeros(length(derivates),1));
for i=1:length(derivates)
    if derivates(i) > 0
        lowerp(i) = p(i).inf;
        upperp(i) = p(i).sup;
    elseif derivates(i) < 0
        lowerp(i) = p(i).sup;
        upperp(i) = p(i).inf;
    elseif derivates(i) == 0
        lowerp(i) = mid(p(i));
        upperp(i) = mid(p(i));
    else
        lowerp(i) = p(i);
        upperp(i) = p(i);
    end
end

% We shrink some parameters. Now we deploy some method to solve our new 
% systems. Place for improvements here. Like if we shrink all parameters,
% we can compute bounds directly.

lowerx = ilspenchbr(A,b,lowep);
upperx = ilspenchbr(A,b,upperp);

iv = infsup(lowerx.inf,upperx.sup);
end