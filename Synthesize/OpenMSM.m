function [freq, A, start, finish] = OpenMSM( filePath )
%OPENBINARY Summary of this function goes here
%   Detailed explanation goes here

    fileid = fopen(filePath, 'r');
    freq = fread(fileid, 1, 'double'); % 8
    size1 = fread(fileid, 1, 'uint32'); % 4
    size2 = fread(fileid, 1, 'uint32'); % 4
    start = fread(fileid, 1, 'uint32'); % 4
    finish = fread(fileid, 1, 'uint32'); % 4
    fread(fileid, [1,5], 'int64'); % 5 * 8bytes for later use.
    A = fread(fileid, [size1, size2], 'double');
    fclose(fileid);
    
end
