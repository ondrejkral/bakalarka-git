function refcheck(solution, refinement, option)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if (in(refinement,solution) | (refinement == solution))
    %disp('REFCHCK: Success!')
else
    %disp('REFCHCK: Fail!')
    disp(option)
    d = inf(solution) - inf(refinement);
    davg = 0; dcount = 0;
    for i = 1:length(d)
        if d(i) > 0
            davg = davg + d(i);
            dcount = dcount + 1;
        end
    end
    
    u = sup(solution) - sup(refinement);
    uavg = 0; ucount = 0;
    for i = 1:length(u)
        if u(i) > 0
            uavg = uavg + u(i);
            ucount = ucount + 1;
        end
    end
    disp('Inf:');
    davg/dcount
    dcount
    disp('Sup');
    uavg/ucount
    ucount
end

end

