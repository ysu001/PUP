function []=findtroisize(valimg, maskimg, N)

[infoval, dataval]=read4dfp(valimg);
[infomask, datamask]=read4dfp(maskimg);

data=dataval(datamask>0);
data=sort(data(:),1,'descend');
t=data(N+1);
fid=fopen('findtroisize_threshold.txt', 'w');
fprintf(fid,'%f\n',t);
fclose(fid);

