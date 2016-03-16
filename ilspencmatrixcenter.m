function V = ilspencmatrixcenter( A, iparameterVector)
%ILSPENCMATRCENTER Compute verified matrix at center of parameter vector.

% Computing center of parameter vector.
parameterCenter = ilspenccenter(iparameterVector);

% Allocation of zero matrix with correct size. 
[m,n] = ilspencmatrixdim(A);
V = intval(zeros(m,n));

% Computing matrix at given center of parameter vector.
for i = 1:length(iparameterVector)
    V = V + ilspencgetak(A,i)*parameterCenter(i);
end
end
