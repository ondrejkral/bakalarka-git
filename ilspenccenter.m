function v = ilspenccenter( ivector )
%ILSPENCCENTER Compute verified center of given interval vector.

v = (inf(ivector) + intval(sup(ivector)))/2;
end

