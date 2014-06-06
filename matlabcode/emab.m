function t=emab(I,NI, m1, s1, m2, s2)


w0=0.001;
w1=0.6;
w2=0.399;
N=sum(NI);
%disp([m1, s1, m2, s2, w0, w1, w2]);
%while (dm1>.1 || dm2>.1 || ds1>.1 || ds2>.1) 
iter=0;
while iter<80,
    p0=1/I(end);
    p1=exp(-1/2*((I-m1)/s1).^2)/sqrt(2*pi*s1^2);
    p2=exp(-1/2*((I-m2)/s2).^2)/sqrt(2*pi*s2^2);
    
    p0old=w0*p0./(w0*p0+w1*p1+w2*p2);
    p1old=w1*p1./(w0*p0+w1*p1+w2*p2);
    p2old=w2*p2./(w0*p0+w1*p1+w2*p2);
    
    m1new=sum(p1old.*I.*NI)/sum(p1old.*NI);
    m2new=sum(p2old.*I.*NI)/sum(p2old.*NI);
    s1new=sqrt(sum(p1old.*(I-m1new).^2.*NI)/sum(p1old.*NI));
    s2new=sqrt(sum(p2old.*(I-m2new).^2.*NI)/sum(p2old.*NI));
    w0new=sum(NI.*p0old)/N;
    w1new=sum(NI.*p1old)/N;
    w2new=sum(NI.*p2old)/N;
    
    m1=m1new;
    m2=m2new;
    s1=s1new;
    s2=s2new;
    w0=w0new;
    w1=w1new;
    w2=w2new;
    
    %disp([m1, s1, m2, s2, w0, w1, w2]);
    iter=iter+1; 
end
    p1=w1*exp(-1/2*((I-m1)/s1).^2)/sqrt(2*pi*s1^2);
    p2=w2*exp(-1/2*((I-m2)/s2).^2)/sqrt(2*pi*s2^2);
    diffp=p2-p1;
    t=find(diffp>0,1);
    if isempty(t), t=I(end); end
