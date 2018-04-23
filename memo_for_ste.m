

%create two sinusoidals signals 
t=1000;
omega=20*pi;

x=zeros(1,t);
y=zeros(1,t);
for tt=1:t
    x(tt)=sin(omega*tt);
end

%create a delay version of signal x
delay=10;
y=[x(delay+1:end) x(1:delay)];

%%%%%%%%%%%%% SIMPLE SYMBOLIZATION SCHEME BASED ON MEAN VALUE %%%%%%%%%
%%%%%% FIRST SIGNAL %%%%%%%
mean1=mean(x);

sx=x;
sx(x > mean1)=1;
sx(x < mean1)=0;


%%%%%% SECOND SIGNAL
mean2=mean(y);

sy=y;
sy(y > mean2)=1;
sy(y < mean2)=0;

%%%%%%%%%%%% run ste %%%%%%%%%%
dif=symbolic_d1transfer_entropy(sx,sy);



%%%%%%%%%%%%% SECOND SYMBOLIZATION SCHEME BASED ON ORDINAL PATTERNS %%%%%%%%%
%%%%%% FIRST SIGNAL %%%%%%%


 L=3;  % # of ordinal patterns
X1=trajectory_matrix(x,L);
X2=trajectory_matrix(y,L); 

 ordinal_patterns=perms([1:L]) ; % # of possible combinations of ordinal patterns equal to L!
 tt=sum([1:L].^2);
 
[ww,ord_pat1]=sort(X1');
ord_pat1=ord_pat1';
[ww,ord_pat2]=sort(X2'); 
ord_pat2=ord_pat2';

% find the ordinal patterns
[row op1]=find((ord_pat1*ordinal_patterns'==tt)==1); 
[row op2]=find((ord_pat2*ordinal_patterns'==tt)==1); 

%%%%%%%%%%%% run ste %%%%%%%%%%
dif=symbolic_d1transfer_entropy(op1,op2);




