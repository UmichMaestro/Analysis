urls = textread('iowa_urls.txt', '%s');
outPath = '/Volumes/WD_Ultra/';
for idx = 1:length(urls)
    url = urls{idx};
    CrackIowa(url, outPath);
end