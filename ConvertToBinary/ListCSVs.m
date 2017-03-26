path = './CSVs'
dirResult = dir(path);
for idx = 1:numel(dirResult)
    item = dirResult(idx);
    fileName = item.name;
    [~,f,e] = fileparts(fileName);
    if strcmp(e,'.csv')
        k = strfind(f, '-');
        freqStr = f(k(end)+1:end);
        freq = str2num(freqStr);
        disp(freq)
    end
end
