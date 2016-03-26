function [ m, n ] = ilspencmatrixdim( A )
%ILSPENCMATRIXDIM Returns matrix dimensions.

global dataModel;
switch(dataModel)
    case '3D'
        m = length(A(:,1,1));
        n = length(A(1,:,1));
    case 'cell'
        m = A{1}(1);
        n = A{1}(2);
end

end

