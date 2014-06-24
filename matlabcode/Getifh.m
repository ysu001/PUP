function [info]=Getifh(filename)
%
% Program to read 4dfp image header. 
%
% Yi Su, Jan 2011

info=[];
%Processing filename
[path,name1,ext1]=fileparts(filename);
[dum,name2,ext2]=fileparts(name1);

if (isempty(ext1)) 
    name=[name1 '.4dfp.ifh'];
elseif (isempty(ext2)||~strcmp(ext2,'.4dfp')||~strcmp(ext1,'.ifh'))
    error('Invalid File');
else
    name=[name2 '.4dfp.ifh'];
end
endian_flag=0;
center_flag=0;
mmppix_flag=0;
% Read the ifh header into memory
fid=fopen(fullfile(path,name),'r');
if (fid==-1), error(['File open error:' fullfile(path,name)]); end
while ~feof(fid)
    line = fgetl(fid);
    if (findstr(line, 'version of keys'))
        [token, remain] = strtok(line,':='); 
        info.version_of_keys=strtrim(remain(3:end));
    elseif (findstr(line, 'conversion program'))
        [token, remain] = strtok(line,':='); 
        info.conversion_program=strtrim(remain(3:end));
    elseif (findstr(line, 'name of data file'))
        [token, remain] = strtok(line,':='); 
        info.name_of_data_file=strtrim(remain(3:end));
    elseif (findstr(line, 'number format'))
        [token, remain] = strtok(line,':='); 
        info.number_format=strtrim(remain(3:end));
    elseif (findstr(line, 'imagedata byte order'))
        [token, remain] = strtok(line,':='); 
        info.imagedata_byte_order=strtrim(remain(3:end));
        endian_flag = 1;
    elseif (findstr(line, 'number of bytes per pixel'))
        [token, remain] = strtok(line,':='); 
        info.number_of_bytes_per_pixel=str2double(remain(3:end));
    elseif (findstr(line, 'number of dimensions'))
        [token, remain] = strtok(line,':='); 
        info.number_of_dimensions=str2double(remain(3:end));
    elseif (findstr(line, 'matrix size [1]'))
        [token, remain] = strtok(line,':='); 
        info.matrix_size(1)=str2double(remain(3:end));
    elseif (findstr(line, 'matrix size [2]'))
        [token, remain] = strtok(line,':='); 
        info.matrix_size(2)=str2double(remain(3:end));
    elseif (findstr(line, 'matrix size [3]'))
        [token, remain] = strtok(line,':='); 
        info.matrix_size(3)=str2double(remain(3:end));
    elseif (findstr(line, 'matrix size [4]'))
        [token, remain] = strtok(line,':='); 
        info.matrix_size(4)=str2double(remain(3:end));
    elseif (findstr(line, 'orientation'))
        [token, remain] = strtok(line,':='); 
        info.orientation=str2double(remain(3:end));
    elseif (findstr(line, 'scaling factor (mm/pixel) [1]'));
        [token, remain] = strtok(line,':='); 
        info.scaling_factor(1)=str2double(remain(3:end));
    elseif (findstr(line, 'scaling factor (mm/pixel) [2]'));
        [token, remain] = strtok(line,':='); 
        info.scaling_factor(2)=str2double(remain(3:end));
    elseif (findstr(line, 'scaling factor (mm/pixel) [3]'));
        [token, remain] = strtok(line,':='); 
        info.scaling_factor(3)=str2double(remain(3:end));
    elseif (findstr(line, 'center'))
        [token, remain] = strtok(line,':='); 
        info.center=sscanf(remain, ':= %f %f %f');
        center_flag=1;
    elseif (findstr(line, 'mmppix'))
        [token, remain] = strtok(line,':='); 
        info.mmppix=sscanf(remain, ':= %f %f %f');  
        mmppix_flag=1;
    end
end
fclose(fid);
if ~endian_flag
    info.imagedata_byte_order='bigendian';
end
if ~mmppix_flag
    info.mmppix=info.scaling_factor.*[1 -1 -1];
end
if ~center_flag
    info.center(1)=info.mmppix(1)*(info.matrix_size(1)-floor(info.matrix_size(1)/2));
    info.center(2)=info.mmppix(2)*(1+floor(info.matrix_size(2)/2));
    info.center(3)=info.mmppix(3)*(1+floor(info.matrix_size(3)/2));
end
    
    