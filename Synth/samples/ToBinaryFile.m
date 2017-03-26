function ToBinaryFile( mat )
%TOBINARYFILE Summary of this function goes here
%   Detailed explanation goes here
    fileid = fopen('out.bin', 'w')
    fwrite(fileid, size(mat,1))
    fwrite(fileid, size(mat,2))
    fwrite(fileid, mat)
    fclose(fileid)

end

