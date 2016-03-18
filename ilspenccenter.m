function v = ilspenccenter( ivector )
%ILSPENCCENTER Compute verified center of given interval vector.

v = (inf(ivector) + sup(ivector))/2;
end

