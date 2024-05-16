function Output= horseshoemcmc(Input)

% n=100;
% K=3;
% 
% X=[ones(n,1)  normrnd(0,1,[n K-1])];
% y=3+2*X(:,2)-2*X(:,3)+normrnd(0,1,[n 1]);
X=Input.zin;
y=Input.yt;
Xout=Input.zout;
delta=Input.delta;
first=Input.first;
[n,K]=size(X);


N=10000;
R=1000;

rue=0;
delta=0.0000;
% Prior distributions
b0=zeros(0,K);
alpha0=0.000;
delta0=0.000;
B=1;
A=1;



%%% save stuff
beta=zeros(K,R+N);
sigma2=zeros(R+N,1);
lambda=zeros(R+N,K);
v=zeros(R+N,K);
tau=zeros(R+N,1);
xi=zeros(R+N,1);

% Initialize
B0=0.001*speye(K); %use a "small variance" (e.g. 0.0001) -> informativ prior 
B0(1:3,1:3)=speye(3);
invB0=inv(B0);
tau(1,:)=1;
lambda(1,:)=ones(1,K);
sigma2(1,1)=1;
alpha1=alpha0+n/2;
b1=zeros(K,1);
%%gibbs sampler
tic
for i=2:R+N

%   if mod(i,100) == 0
%          disp([num2str(i) ' Simulations'])
%          disp([num2str( sum(ind)) 'active'])       
%   end
%   
  
  delta1=delta0+(y-X*beta(:,i-1))'*(y-X*beta(:,i-1))/2;  
  sigma2(i,1)=1/gamrnd(alpha1,1/delta1,1);
  
%   B1=inv((1/sigma2(i,1))*X'*X+ invB0);
%   b1=B1*((1/sigma2(i,1))*X'*y);%+solve(B0)*b0
%   B1chol=chol(B1,'lower');  
%   beta(:,i)= b1 +B1chol*randn(K,1);

  


if rue==1
    
  B1=((1/sigma2(i,1))*X'*X+ invB0);
  
  b1=B1\((1/sigma2(i,1))*X'*y);%+solve(B0)*b0
 

  B1chol=chol(B1,'lower');
  beta(:,i)= b1 +B1chol'\randn(K,1); %(range{ii},i)

else
    
  Xtilde=(1/sqrt(  sigma2(i,1)))*X;
  ytilde=y*(1/sqrt(  sigma2(i,1)));
  V=  diag(B0.^.5).*randn(K,1);
  q=randn(n,1);
  w=Xtilde*V+q;
  ind=diag(B0)>delta;
  Xtildeadj=Xtilde(:,ind);
  B1=(speye(n)+Xtildeadj*B0(ind,ind)*Xtildeadj');
  u=B1\(ytilde-w);   % u=inv(B1)*(ytilde-w)
  beta(:,i)=B0(:,ind)*Xtildeadj'*u+V;%(range{ii},i)
end

 for j = 1:K
    v(i,j)=1/gamrnd(1,1/(1/B^2+1/lambda(i-1,j)),1);%v[i,]
    lambda(i,j)=1/gamrnd(1,1/(1/v(i,j)+(beta(j,i)^2)/(2*tau(i-1,1)) ),1) ;% lambda[i,]
 end

  xi(i,1)=1/gamrnd(1,1/(1/A^2+1/tau(i-1,1) ),1);

  tau(i,1)=1/gamrnd((K+1)/2,1/( 1/xi(i,1)+0.5*sum(beta(:,i)'.^2 ./lambda(i,:)) ),1);
  invB0(1:(K+1):end)= ((1./(lambda(i,:)*tau(i,1))));
  B0(1:(K+1):end)=(((lambda(i,:)*tau(i,1))));

  
  
end

%delete burn-in period
beta = beta(:,(R+1):(R+N));

betamean=mean(beta,2);

Output.betamean=betamean;
