function V = ilspencmakeb( b, position, p)
%ILSPENCMAKEB Returns vector in current representation and with par = 1.
%   Depends on which representation is switched (3D, cell, ...). Returns
%   vector on parameter position given. However, values are intervals.
global dataModel;
switch(dataModel)
    case 'cell'
        V{1} = [length(b); 1];
        V{2} = intval(zeros(2,length(b)));
        for i = 1:length(b)
            V{2}(:,i) = [intval(i), b(i)];
        end
        for i = 2:length(p)
            V{i+1} = [1, 0];
        end
                
    case '3D'
        V = intval(zeros(length(b),length(p)));
        V(:,1) = b;
    otherwise
        % This case shouldn't occurs beacuse this exception is handled on
        % higher layer.
end

end

