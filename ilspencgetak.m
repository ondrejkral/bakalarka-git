function V = ilspencgetak( A, k )
%ILSPENCGETAK Returning matrix of coeficients for k-th parameter.
%   Additional layer for testing repesentation of parameters' linear 
%   dependencies. Because of MATLAB's lazy copying this shouldn't be
%   significant performance overhead.

global dataModel;
switch(dataModel)
    case '3D'
        V = intval(A(:,:,k));
    case 'cell'
        % obtaining matrix dimension and other meta-data
        m = A{1}(1);
        n = A{1}(2);
        
        V = zeros(m,n);
        Ak = A{k+1};
        for i = 1:length(Ak(1,:))
            column = Ak(:,i); 
            V(column(1),column(2)) = column(3);
        end
        
        V = intval(V);
        
end
end

