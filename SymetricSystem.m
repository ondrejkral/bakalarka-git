function [A,b,p] = SymetricSystem(coefMagnitudeMult, dimension, R)
parametersCount = (dimension*dimension - dimension)/2 + dimension;
randomACenter = -10 + 20*rand(1,parametersCount);
randomACenter = 10*coefMagnitudeMult + randomACenter;
p = midrad(randomACenter,R);
p(1) = intval('1');
randomBCenter = -10 + 20*rand(1,dimension);

global dataModel;

switch(dataModel)
    case 'cell'
        A{1} = [dimension; dimension; 0;];
        b{1} = [dimension; 0];
        for i = 1:dimension
            b{2}(:,i)= [i; randomBCenter(i)];
        end
        
        for j = 1:(parametersCount - 1)
            b{2 + j} = [1; 0];
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
        b = zeros(dimension,parametersCount);
        b(:,1) = randomBCenter;
        
        A = zeros(dimension, dimension, parametersCount);
        
        parameterIndex = 1;
        for j = 1:dimension
            for i = 1:j
                if (i == j)
                    A(i,i,parameterIndex) = 1; 
                else
                    A(i,j,parameterIndex) = 1;
                    A(j,i,parameterIndex) = 1;
                end
                parameterIndex = parameterIndex + 1;
            end
        end
    otherwise
        %invalid model
end
end
