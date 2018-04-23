function dif=symbolic_dtransfer_entropy(x,y,delay)

%This m - file computes the symbolic tranfer entropy between tw symbolic
%time series
%Reference:
%Staniek & Lehnertz,"Symbolic trasnfer entropy", PHYSICAL REVIEW LETTERS,
%2008

%INPUT : symbolix time series x,y 
%        delay -> delay parameter

%OUTPUT: dif = tentxy - tentyx

%DIMITRIADIS STAVROS  10/2012
%Dominik Storhas  7/2013 (correction of pxy estimation - line 153)
%Technical University Munich (Germany) and 
%Macquarie University Sydney (Australia)


s=1;
uni=unique([x,y]);
nosym=length(uni);


for k=1:length(x)
    r=find(x(k)==uni);
    x(k)=r;
    
    r=find(y(k)==uni);
    y(k)=r;
end


%ESTIMATE TRANSFER ENTROPY Y -> X
dif=zeros(1,delay);

for d=1:delay
    
%estimate pxxy
pxxy=zeros(nosym,nosym,nosym);

for k=d:length(x)-d
    pxxy(x(k+s),x(k),y(k+s-d))=pxxy(x(k+s),x(k),y(k+s-d))+1;
end

sum1=sum(sum(sum(pxxy)));
pxxy=pxxy/sum1;

%estimate pxx
    pxx=zeros(nosym,nosym);

    for k=1:length(x)-s
        pxx(x(k+s),x(k))=pxx(x(k+s),x(k))+1;
    end


sum1=sum(sum(pxx));
pxx=pxx/sum1;



%estimate pxy
    pxy=zeros(nosym,nosym);

    for k=1:length(x)
        pxy(x(k),y(k))=pxy(x(k),y(k))+1;
    end

sum1=sum(sum(pxy));
pxy=pxy/sum1;

%estimate px
px=zeros(1,nosym);


for k=1:length(x)
    px(x(k))=px(x(k))+1;
end

sum1=sum(px);
px=px/sum1;


%TRANSFER ENTROPY Y -> X
tentyx=0;

count=0;
for k=1:nosym
    for l=1:nosym
        for m=1:nosym
            count=count + 1;
            num=pxxy(k,l,m)*px(l);
            dem=pxy(l,m)*pxx(k,l);
            tentyx(count)=pxxy(k,l,m)*log2(num/dem);
        end
    end
end

%eliminate NaNs 
tentyx=sum(tentyx(find(isnan(tentyx)==0)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ESTIMATE TRANSFER ENTROPY X -> Y


%estimate pyyx
pyyx=zeros(nosym,nosym,nosym);

for k=d:length(y)-d
    pyyx(y(k+s),y(k),x(k+s-d))=pyyx(y(k+s),y(k),x(k+s-d))+1;
end

sum1=sum(sum(sum(pyyx)));
pyyx=pyyx/sum1;

%estimate pyy
pyy=zeros(nosym,nosym);

for k=1:length(y)-s
    pyy(y(k+s),y(k))=pyy(y(k+s),y(k))+1;
end


sum1=sum(sum(pyy));
pyy=pyy/sum1;


%estimate py
py=zeros(1,nosym);


for k=1:length(y)
    py(y(k))=py(y(k))+1;
end

sum1=sum(py);
py=py/sum1;




%TRANSFER ENTROPY X -> Y
tentxy=0;

count=0;
for k=1:nosym
    for l=1:nosym
        for m=1:nosym
            count=count + 1;
            num=pyyx(k,l,m)*py(l);
            dem=pxy(m,l)*pyy(k,l);
            tentxy(count)=pyyx(k,l,m)*log2(num/dem);
        end
    end
end

%eliminate NaNs 
tentxy=sum(tentxy(find(isnan(tentxy)==0)));


dif(d)=tentxy - tentyx;

end

plot(1:1:delay,dif)
xlabel('delay in ms')
ylabel('delay tranfer entropy')


