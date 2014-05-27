function [status]=emabv1(TOFfn, m1, s1, m2, s2)
% Program for TOF-MRA segmentation to obtain air and bone mask.
% USAGE:
%   [status]=emabv1(TOFfn, m1, s1, m2, s2);
%   Empirical values for inputs:
%   m1=6; s1=4; m2=100; s2=50;
%
%   Algorithm adapted from Wilson & Noble, IEEE TMI 1999 18(10) 938-945. 
%   Yi Su, 11/11/2011
%

[infoTOF,dataTOF]=read4dfp(TOFfn);
sz=size(dataTOF);
ns=sz(3);
dataT=dataTOF;
for i=1:ns
    tmp=dataTOF(:,:,i);
    I=(0:max(tmp(:)))';
    pall=histc(tmp(:),I);
    T(i)=emab(I,pall,m1,s1,m2,s2);
    tmp=dataT(:,:,i);
    tmp1=tmp;
    tmp1(tmp>T(i))=0;
    tmp1(tmp<=T(i))=1;
    dataT(:,:,i)=tmp1;
end
infoT=infoTOF;
infoT.conversion_program='matlab';
infoT.name_of_data_file='TOF_AB';
status=write4dfp('TOF_AB',infoT,dataT);
