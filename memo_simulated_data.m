
N=2^12; 
P=3; 


myswitch=5;%if positive direction is from channel 1 to channel 2;

t=0:0.05:1;

ste=zeros(1,length(t));

nn=20; % # of realizations
dif_ab=zeros(nn,length(t));

matlabpool; % for multi-cores

tic
for k1=1:nn
    parfor k2=1:length(t)
    [k1 k2]

data=simuldata(N,P,t(k2),myswitch);

signal1=data(:,1);
signal2=data(:,2);


[b,a]=butter(2,0.1,'low'); 
filt1=filter(b,a,signal1');
filt2=filter(b,a,signal2');

noprot=4;
iterations=10;

m=3;
tau=10;

[esignal1] = embeddelay(signal1, m, tau);
[esignal2] = embeddelay(signal2, m, tau);

tic,[prot,class]=Vector_Quantization([esignal1;esignal2],noprot,iterations);,toc


len=length(class)/2;
ct1=class(1:len);
ct2=class(len+1:end);

%dif=symbolic_d1transfer_entropy(ct1,ct2);
nodel=100;
dif=symbolic_dtransfer_entropy(ct1,ct2,nodel);

%%% find the maximum value across the ones with zscore > 2
 % dif_val=max(dif(find(abs(zscore(dif)) > 2)));
 
  
  %%%%%%% surrogate analysis %%%%%%%%%%%%%%%
  nosur=1000;
  surrogates=zeros(nodel,nosur);
  
  for kk=1:nosur
      rr=randperm(length(ct2));
      shuffle2=ct2(rr);
      dif_sur=symbolic_dtransfer_entropy(ct1,shuffle2,nodel);
      surrogates(:,kk)=dif_sur;
  end
  
  pval=zeros(1,nodel);
  for kk=1:nodel
    pval(kk)=(length(find(surrogates(kk,:) > dif(kk))) + 1 )/nosur;
  end

  
  [val delay]=min(pval);
  
   if val < 0.05
      rr=find(pval==val); %% find dste values with equal min pval
      
      if isempty(rr)==1
          dif_ab(k1,k2)=dif(delay);
      elseif isempty(rr)==0
          no_pos=0;
          no_neg=0;
          no_pos=length(find(dif(rr) > 0));
          no_neg=length(find(dif(rr) < 0));
          
          if no_pos > no_neg
              dif_ab(k1,k2)=max(dif(rr));
          elseif no_neg > no_pos
              dif_ab(k1,k2)=min(dif(rr));
          end
      end
  end
 

    end
end

toc

%%% quantify the accuracy (%) of correct detection %%%%%%%%
[r c]=find(isinf(dif_ab)==1);

for k=1:length(r)
    dif_ab(r(k),c(k))=0;
end

perfomance=sum(dif_ab > 0)./nn;


%%% get the mean and std of positive values

pos=(dif_ab.*(dif_ab > 0));
total=sum(dif_ab > 0);
meanval=sum(pos)./total;
[r c]=find(dif_ab);


for k=1:length(t)
    stdval(k)=std(nonzeros(pos(:,k)));
end

figure(1);errorbar(meanval,stdval)
figure(2);plot(performance)



