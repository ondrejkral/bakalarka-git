function [A,b,p] = SymetricSystem(coefMagnitudeMult, dimension, R)
parametersCount = (dimension*dimension - dimension)/2 + dimension;
randomACenter = -10 + 20*rand(1,parametersCount);
randomACenter = 10*coefMagnitudeMult + randomACenter;
p = midrad(randomACenter,R);
randomBCenter = -10 + 20*rand(1,dimension);

global dataModel;

switch(dataModel)
    case 'cell'
        A{1} = [dimension; dimension; 0;];
        b{1} = [dimension; 0];
        for i = 1:dimension
            b{2}(:,i)= [i; randomBCenter(i)];
        end
        
        for j = 3:parametersCount
            b{j} = [1; 0];
        end
            
        parameterIndex = 1;
        for j = 1:dimension
            for i = 1:j
                if (i == j)
                    A{parameterIndex+1} = [i; i; 1;]; 
                else
                    A{parameterIndex+1} = [i, j; j, i; 1, 1;];
                end
                parameterIndex = parameterIndex + 1;
            end
        end
    case '3D'
        %%%
    otherwise
        %invalid model
end
end
