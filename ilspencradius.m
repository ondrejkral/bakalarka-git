function v = ilspencradius( ivector )
%ILSPENCRADIUS Compute verified radius of given interval vector.

v = (sup(ivector) - intval(inf(ivector)))/2;
end

