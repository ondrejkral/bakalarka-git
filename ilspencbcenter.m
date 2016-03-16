function v = ilspencbcenter( b, iparameterVector )
%ILSPENCBCENTER Compute verified right side of linear system at center of 
% parametric vector.

% Computing center of parameter vector.
parameterCenter = ilspenccenter(iparameterVector);

% Allocation of zero vector with correct size. 
l = ilspencvectordim(b);
v = intval(zeros(l,1));

% Computing vector at given center of parameter vector.
for i = 1:length(iparameterVector);
    v = v + ilspencgetbk(b,i)*parameterCenter(i);
end
end

