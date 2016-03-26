function [A,b,p] = SymetricSystem(coefMagnitudeMult, dimension, R)
% Number of parameters needed.
parametersCount = (dimension*dimension - dimension)/2 + dimension;

% Random parameters center.
random = -10 + 20*rand(parametersCount,1);
% randomACenter = 10*coefMagnitudeMult + randomACenter;
% parameters vector
p = midrad(random,R);

% random vector
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
                    p(parameterIndex) = p(parameterIndex) + 10*dimension;
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
                     p(parameterIndex) = p(parameterIndex) + 10*dimension;
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
