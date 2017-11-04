path = '.'; % set path
dirResult = dir(path);
for idx = 1:numel(dirResult)
    item = dirResult(idx);
    fileName = item.name;
    [~,f,e] = fileparts(fileName);
    if strcmp(e,'.csv')
        k = strfind(f, '-');
        freqStr = f(k(end)+1:end);
        freq = str2num(freqStr);
        freq = freq(1)
        
        
%         fileid = fopen('out.bin', 'w')
%     fwrite(fileid, size(mat,1), 'uint')
%     fwrite(fileid, size(mat,2), 'uint')
%     fwrite(fileid, mat, 'double')
%     fclose(fileid)
    end
end
