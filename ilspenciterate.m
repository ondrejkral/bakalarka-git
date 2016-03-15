function v = ilspenciterate( A, b, p, iterations, option )
%ILSPENCITERATE Divide parameters' space into several sub-spaces and
% solve them recursively.

% bottom of recursion
if iterations == 0
    switch option
        case 'BS'
            v = ilspencbauerskeel(A,b,p);
        case 'HBR'
            v = ilspenchbr(A,b,p);
        case 'SKALNA'
            v = ilspencresidual(A,b,p, 'SKALNA');
        case 'RUMP'
            v = ilspencresidual(A,b,p, 'RUMP');
        otherwise
            % invalid option
            v = NaN;
    end
else

    % find largest interval
    max = 0; maxindex = 0; pi = intval(0);
    for i = 1:length(p)
        if (sup(p(i)) - inf(p(i))) > max
            max = sup(p(i)) - inf(p(i));
            maxindex = i;
            pi = p(i);
        end
    end

    % split largest interval
    p1 = p; p1(maxindex) = infsup(inf(pi), inf(pi) + max/2);
    p2 = p; p2(maxindex) = infsup(inf(pi) + max/2, sup(pi));
    
    % recursion
    disp(iterations);
    v = hull(ilspenciterate(A, b, p1, iterations - 1, option),ilspenciterate(A, b, p2, iterations - 1, option));
end
end

