function [ A, b, p ] = toeplitzsystem( coefmagnitude, radius, matrixmagnitude)

global dataModel;

% random center of parameters
random = -coefmagnitude + (2*coefmagnitude)*rand(2*matrixmagnitude - 1,1);
p = intval(zeros(1,2*matrixmagnitude - 1));
% first parameter is set as
p(1) = intval('1');
for i = 2:(2*matrixmagnitude - 1)
    p(i) = midrad(random(i),radius);
end
% ///

switch(dataModel)
    case '3D'
        A = zeros(matrixmagnitude,matrixmagnitude,2*matrixmagnitude - 1);
        
        c = zeros(matrixmagnitude,1); c(1) = 1;
        r = zeros(matrixmagnitude,1);

        for i = (matrixmagnitude+1):(2*matrixmagnitude - 1)
            c = circshift(c,1);
            A(:,:,i) = toeplitz(c,r);
        end

        c = zeros(matrixmagnitude,1); c(1) = 1;
        r = zeros(matrixmagnitude,1); r(1) = 1;
        for i = 1:matrixmagnitude
            A(:,:,i) = toeplitz(c,r);
            c(1) = 0;
            r = circshift(r,1);
        end

        b = zeros(matrixmagnitude,2*matrixmagnitude-1);
        b(:,1) = -coefmagnitude + (2*coefmagnitude)*rand(matrixmagnitude,1);
    case 'cell'
        A{1} = [ matrixmagnitude; matrixmagnitude; 0;];
        b{1} = [matrixmagnitude; 0;];
        
        for i = 1:matrixmagnitude
            A{i+1} = zeros(3,matrixmagnitude - i + 1);
            for j = 1:(matrixmagnitude - i + 1)
                A{i+1}(:,j)=[j; (j+i-1); 1;];
            end
        end
        for i = 2:matrixmagnitude
            A{matrixmagnitude + i} = zeros(3,matrixmagnitude - i + 1);
            for j = 1:(matrixmagnitude - i + 1)
                A{matrixmagnitude + i}(:,j)=[(j+i-1); j; 1;];
            end
        end
        randomvector = -coefmagnitude + (2*coefmagnitude)*rand(matrixmagnitude,1);
        for j = 1:matrixmagnitude 
            b{2}(:,j) =[ j; randomvector(j);] ;
        end
        for i = 3:(2*matrixmagnitude+1)
            b{i} = [1; 0;];
        end
            
end                
end