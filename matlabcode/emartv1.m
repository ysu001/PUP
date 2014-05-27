function [status]=emartv1(TOFfn, m1, s1, m2, s2, m0, s0)
% Program for TOF-MRA segmentation to obtain arterial mask.
% USAGE:
%   [status]=emartv1(TOFfn, m1, s1, m2, s2, m0, s0);
%   Empirical values for inputs:
%   m1=6; s1=4; m2=100; s2=50; m0=300; s0=100;
%
%   Algorithm adapted from Wilson & Noble, IEEE TMI 1999 18(10) 938-945. 
%   Yi Su, 10/28/2011
%

[infoTOF,dataTOF]=read4dfp(TOFfn);
sz=size(dataTOF);
ns=sz(3);
dataT=dataTOF;
for i=1:ns
    tmp=dataTOF(:,:,i);
    I=(0:max(tmp(:)))';
    pall=histc(tmp(:),I);
    T1(i)=emart(I,pall,m1,s1,m2,s2);
    T2(i)=emart2(I,pall,m1,s1,m2,s2,m0,s0);
    T(i)=max(T1(i),T2(i));
    tmp=dataT(:,:,i);
    tmp(tmp<T(i))=0;
    dataT(:,:,i)=tmp;
end
infoT=infoTOF;
infoT.conversion_program='matlab';
infoT.name_of_data_file='TOF_ART';
status=write4dfp('TOF_ART',infoT,dataT);
