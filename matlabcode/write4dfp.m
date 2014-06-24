function [status]=write4dfp(filename, info, data)

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
    error(['Invalid File Name' fullfile(path,name)]);
elseif (strcmp(ext1,'.img')&&strcmp(ext2,'.4dfp'))
    name=[name2 '.4dfp.img'];
    nameroot=name2;
elseif (strcmp(ext1,'.ifh')&&strcmp(ext2,'.4dfp'))
    name=[name2 '.4dfp.img'];
    nameroot=name2;
else
    error(['Invalid File Name: ' fullfile(path,name)]);
end

%Write .4dfp.img file
if strcmp(info.imagedata_byte_order,'bigendian')
    BO='b';
else
    BO='l';
end
fid=fopen(fullfile(path,name),'w',BO);
if (fid==-1), error(['File open error:' fullfile(path,name)]); end
if info.matrix_size(4)>1
    data=permute(data,[2 1 3 4]);
else
    data=permute(data,[2 1 3]);
end
status=fwrite(fid,data,'float32',BO);
fclose(fid);

%Write .4dfp.ifh file
Writeifh(fullfile(path,nameroot),info);

%Write .4dfp.hdr file
system(['ifh2hdr ' fullfile(path,nameroot)]);
return;