function []=Writeifh(filename, ifhdr)

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

% Write the ifh header file
fid = fopen(fullfile(path,name),'w');
if (fid==-1), error(['File open error:' fullfile(path,name)]); end
fprintf(fid, 'INTERFILE\t:=\n');
fprintf(fid, 'version of keys\t:= %s\n', ifhdr.version_of_keys);
fprintf(fid, 'number format\t\t:= %s\n', ifhdr.number_format);
fprintf(fid, 'conversion program\t:= %s\n', ifhdr.conversion_program);
fprintf(fid, 'name of data file\t:= %s\n', ifhdr.name_of_data_file);
fprintf(fid, 'number of bytes per pixel\t:= %d\n', ifhdr.number_of_bytes_per_pixel);
fprintf(fid, 'imagedata byte order\t:= %s\n', ifhdr.imagedata_byte_order);
fprintf(fid, 'orientation\t\t:= %d\n', ifhdr.orientation);
fprintf(fid, 'number of dimensions\t:= %d\n', ifhdr.number_of_dimensions);
for i=1:4
    fprintf(fid, 'matrix size [%d]\t:= %d\n', i, ifhdr.matrix_size(i));
end
for i=1:3
    fprintf(fid, 'scaling factor (mm/pixel) [%d]\t:= %f\n', i, ifhdr.scaling_factor(i));
end
fprintf (fid, 'mmppix\t:= %10.6f%10.6f%10.6f\n', ifhdr.mmppix(1), ifhdr.mmppix(2), ifhdr.mmppix(3));
fprintf (fid, 'center\t:= %10.4f%10.4f%10.4f\n', ifhdr.center(1), ifhdr.center(2), ifhdr.center(3));

fclose(fid);