%$Header: /data/petsun4/data1/solaris/csh_scripts/RCS/alpay.m,v 1.1 2008/08/15 22:48:14 avi Exp $
%$Log: alpay.m,v $
% Revision 1.1  2008/08/15  22:48:14  avi
% Initial revision
%
%Copyright KILICHAN GURLEYIK March 2003
%Modified the original code written by I. Alpay Ozcan with his consent 
%Obtain the consent of both authors before using this software
%It is illegal to distribute this software without the consent of the authors

figure;

%Orientation of the slice 
slplane=0;
while (slplane < 1 | slplane>3)
disp('Please enter a number between 1 and 3 for slice orientation');
slplane=input('Enter the slice plane 1=xy(Trans) 2=xz(Corr) 3=yz(Sag) : ');
end

%the number of slices for each view
maxnoslices=[nz ny nx];
maxnoslices=maxnoslices(slplane);

rem=['Please enter a number between 1 and ' num2str(maxnoslices)];
slno=0;
while (slno < 1 | slno > maxnoslices)
    disp(rem);
    slno=input('Enter the slice number: ');
end

%assign the matrix for the background image


switch slplane
case 1 %xy plane
    sagit=reshape(testall(:,:,slno),[nx ny]);
    rem='xy plane';
    %whisker=bigwhisker(bigwhisker(:,3)==slno,:);
    %evalstr='plot(xp,yp,''b'');plot(xs,ys,''g'');';
  
case 2  %xz plane
    sagit=reshape(testall(:,slno,:),[nx nz]);
    rem='xz plane';
    %whisker=bigwhisker(bigwhisker(:,2)==slno,:);
    %evalstr='plot(xp,zp,''b'');plot(xs,zs,''g'');';
    
case 3  %yz plane
    sagit=reshape(testall(slno,:,:),[ny nz]);
    rem='yz plane';
    %whisker=bigwhisker(bigwhisker(:,1)==slno,:);
    %evalstr='plot(yp,zp,''b'');plot(ys,zs,''g'');';        
    
end
%assign the whisker data
whisker=bigwhisker(bigwhisker(:,4-slplane)==slno,:);

%put the background image
%imagesc(sagit',[min(sagit(:)) max(sagit(:))]);colormap(gray);
imagesc(sagit',[0 600]);colormap(gray);
set(gca,'YDir','normal');
title([rem '  slice ' num2str(slno)]);
axis image;
hold on;

scal=7;

%prepare the whiskers to be plotted
xp=[whisker(:,1)-scal*whisker(:,4)/2 ...
        whisker(:,1)+scal*whisker(:,4)/2 ...
        repmat(NaN,length(whisker),1)]';

yp=[whisker(:,2)-scal*whisker(:,5)/2 ...
        whisker(:,2)+scal*whisker(:,5)/2 ...
        repmat(NaN,length(whisker),1)]';

zp=[whisker(:,3)-scal*whisker(:,6)/2 ...
        whisker(:,3)+scal*whisker(:,6)/2 ...
        repmat(NaN,length(whisker),1)]';

xs=[whisker(:,1)-scal*whisker(:,7)/2 ...
        whisker(:,1)+scal*whisker(:,7)/2 ...
        repmat(NaN,length(whisker),1)]';

ys=[whisker(:,2)-scal*whisker(:,8)/2 ...
        whisker(:,2)+scal*whisker(:,8)/2 ...
        repmat(NaN,length(whisker),1)]';

zs=[whisker(:,3)-scal*whisker(:,9)/2 ...
        whisker(:,3)+scal*whisker(:,9)/2 ...
        repmat(NaN,length(whisker),1)]';

%plot the whiskers
evalstr=['plot(xp,yp,''y'');%plot(xs,ys,''g'');';
    'plot(xp,zp,''y'');%plot(xs,zs,''g'');';
    'plot(yp,zp,''y'');%plot(ys,zs,''g'');'];

eval(evalstr(slplane,:));

disp('To make another whisker plot please type ''alpay''');
%plot3(xs,ys,zs,'g');




