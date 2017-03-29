function WriteSustainPoints( path, start, finish )
%WRITESUSTAINPOINTS Summary of this function goes here
%   Detailed explanation goes here
    fileid = fopen(path, 'r+');
    fseek(fileid, 16, 'bof'); % 8: frequency, 4: duration, 4: partials
    fwrite(fileid, start, 'uint32'); % 4
    fwrite(fileid, finish, 'uint32'); % 4
    fclose(fileid);
end

