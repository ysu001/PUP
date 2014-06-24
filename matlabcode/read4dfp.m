function [info, data]=read4dfp(filename)
%
% Program to read 4dfp format data. 
%
% Yi Su, Jan 2011


%Processing filename
[path,name1,ext1]=fileparts(filename);
[dum,name2,ext2]=fileparts(name1);

if (isempty(ext1)) 
    name=[name1 '.4dfp.img'];
    nameroot=name1;
elseif (isempty(ext2)&&strcmp(ext1,'.4dfp'))
    name=[name1 '.4dfp.img'];
    nameroot=name1;
elseif (isempty(ext2)&&~strcmp(ext1,'.4dfp'))
    error(['Invalid File Name: ' fullfile(path,name)]);
elseif (strcmp(ext1,'.img')&&strcmp(ext2,'.4dfp'))
    name=[name2 '.4dfp.img'];
    nameroot=name2;
elseif (strcmp(ext1,'.ifh')&&strcmp(ext2,'.4dfp'))
    name=[name2 '.4dfp.img'];
    nameroot=name2;
else
    error(['Invalid File Name: ' fullfile(path,name)]);
end

%Read header info
[info]=Getifh(fullfile(path,nameroot));

%Read image data
if strcmp(info.imagedata_byte_order,'bigendian')
    BO='b';
else
    BO='l';
end
fid=fopen(fullfile(path,name),'r',BO);
if (fid==-1), error(['File open error:' fullfile(path,name)]); end
data=single(zeros(info.matrix_size(2), info.matrix_size(1), info.matrix_size(3), info.matrix_size(4)));
for frame=1:info.matrix_size(4)
    [tmp, count]=fread(fid,info.matrix_size(1)*info.matrix_size(2)*info.matrix_size(3),'*float32',0,BO);
    tmp=reshape(tmp,[info.matrix_size(1), info.matrix_size(2), info.matrix_size(3)]);
    tmp=permute(tmp,[2 1 3]);
    data(:,:,:,frame)=tmp;
end
fclose(fid);