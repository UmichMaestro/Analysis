function ToBinaryFile( mat )
%TOBINARYFILE Summary of this function goes here
%   Detailed explanation goes here
    fileid = fopen('out.bin', 'w')
    fwrite(fileid, size(mat,1), 'uint')
    fwrite(fileid, size(mat,2), 'uint')
    fwrite(fileid, mat, 'double')
    fclose(fileid)

end

