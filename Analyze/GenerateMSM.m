function GenerateMSM( path )
dirResult = dir(path);
for idx = 1:numel(dirResult)
    item = dirResult(idx);
    fileName = item.name;
    [~,f,e] = fileparts(fileName);
    if contains(f, 'wav_orig')
        disp(fileName)
        [s,fs] = audioread(sprintf('%s/%s',item.folder,item.name));
        [frq, A] = Analysis({s}, fs);
        outName = sprintf('%s/%s.msm',item.folder,f(1:end-9));
        WriteToFile(frq, A, outName);
        disp(frq)
    end
end

end

function WriteToFile( freq, A, fileName )
%TOBINARYFILE Summary of this function goes here
%   Detailed explanation goes here
    fileid = fopen(fileName, 'w');
    fwrite(fileid, freq, 'double'); % 8
    fwrite(fileid, size(A,1), 'uint32'); % 4
    fwrite(fileid, size(A,2), 'uint32'); % 4
    fwrite(fileid, [0,0,0,0,0,0], 'int64'); % 6 * 8bytes for later use.
    fwrite(fileid, A, 'double');
    fclose(fileid);
end


%  for idx = 1:numel(dirinfo)
%      item = dirinfo(idx);
%      path = sprintf('%s/%s', item.folder, item.name);
%      GenerateBinaryTable(path)
%  end
