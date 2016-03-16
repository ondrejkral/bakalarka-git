function v = ilspencvectordim( b )
%ILSPENCVECTORDIM Return length of vector b

global dataModel;
switch(dataModel)
    case '3D'
        v = length(b(:,1));
    case 'cell'
        v = b{1}(1);
end
end

