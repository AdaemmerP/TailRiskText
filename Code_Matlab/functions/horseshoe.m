function Output= horseshoe(Input)

% n=100;
% K=30;
% 
% X=[ones(n,1)  normrnd(0,1,[n K-1])];
% y=3+2*X(:,2)-2*X(:,3)+normrnd(0,1,[n 1]);
X=Input.zin;
y=Input.yt;
Xout=Input.zout;
delta=Input.delta;
[n,K]=size(X);
quant=Input.quant;
first=Input.first;




% Prior distributions
alpha0=0.00;
delta0=0.00;



% Initialize

betamean=zeros(K,1);
oldbetmean=ones(K,1)*10000;
betasigma=eye(K);
lambda=ones(K-first,1);
v=ones(K-first,1);
tau_beta=0.01;
invB0= eye(K)*10;
invB0(1:first,1:first)=eye(first)*0.0001;
alpha1=alpha0+1.5*n;
tau = 2/(quant*(1-quant));
theta = (1-2*quant)/(quant*(1-quant)); 
mean_inv_vi=ones(n,1);
mean_vi=ones(n,1);
XsigmaX=ones(n,1);
%%
check=0;
count=0;
while check==0
    count=count+1;

   
      
     if delta<1 && delta>0
      
          s_tinv = zeros(n,1);
          at = zeros(n,1);
          bt = zeros(n,1);

    for t = 1:n          
        temp =mean_vi(t)+ (mean_inv_vi(t)*( y(t)-  X(t,:)*betamean)^2+2*theta*(X(t,:)*betamean-y(t))...
         +mean_vi(t)*theta^2+ mean_inv_vi(t)*X(t,:).^2*diag(betasigma(1:end,1:end)) )/(2*tau);
        if t == 1
            at(t,1) = 0 + 1.5;
            bt(t,1) = 0 + temp;
        else
            at(t,1) = delta*at(t-1,:) + 1.5;
            bt(t,1) = delta*bt(t-1,:) + temp;
        end
        s_tinv(t,1) = at(t,1)./bt(t,1);
    end
    % Smooth volatilities
    phi = zeros(n,1); phi(n,1) = at(n,1)./bt(n,1);
    for t=n-1:-1:1
        phi(t,:) = (1-delta)*s_tinv(t,:) + delta*phi(t+1,:);
    end
    %sigma_t = 1./phi;      
    invsigma=diag(phi);
      alpha1= at(n,1);
      delta1=bt(n,1);
            
    elseif delta==0
     invsigma=eye(n); 
     else     
           Sum=zeros(n,1);
     for nn=1:n
   % Sum(nn)= y(nn)^2 .*mean_inv_vi(nn) -mean_inv_vi(nn)*2*y(nn)*X(nn,:)*betamean...
         %+ 2*X(nn,:)*betamean+mean_vi(nn)-2*y(nn)+mean_inv_vi(nn)*X(nn,:).^2*(betamean.^2+diag(betasigma(1:end,1:end)) ) ; 
      Sum(nn)=mean_inv_vi(nn)*( y(nn)-  X(nn,:)*betamean)^2+2*theta*(X(nn,:)*betamean-y(nn))...
         +mean_vi(nn)*theta^2+ mean_inv_vi(nn)*X(nn,:).^2*diag(betasigma(1:end,1:end)) ; 
     end
    Sum(Sum<0)=0;
      delta1=delta0+sum(mean_vi)+0.5*(1/tau)*sum(Sum);
      invsigma=eye(n)*alpha1/delta1;  
    end
    
    
  for nn=1:n
 XsigmaX(nn,1)=X(nn,:)*betasigma*X(nn,:)';
 end
 
 b=( (y-X*betamean).^2+ XsigmaX ).*diag( invsigma)/tau;
 a= 2*diag( invsigma) +diag( invsigma)* theta^2 /tau;
 mean_inv_vi= meanGIG2(.5,a,b,-1);
 mean_vi= meanGIG2(.5,a,b,1);
 
  if count>1
    oldbetmean=  betamean; 
  end
  

   betasigma=( X'*invsigma*diag(mean_inv_vi)*X/tau+ invB0)\eye(K);
 
    betamean=betasigma*(X'*invsigma*diag(mean_inv_vi)*(y-theta./mean_inv_vi)/tau);

  
    for j =1:K-first
    v(j)=1/(1+lambda(j));
    lambda(j)=1/(v(j)+tau_beta*(betamean(j+first)^2+betasigma(j+first,j+first))*0.5);
    end
  xi=1/(1+tau_beta);  
  tau_beta=(K+1-first)/(2*xi+sum((betamean(1+first:end).^2+diag(betasigma(1+first:end,1+first:end))).*lambda));
  
    shrink= lambda*tau_beta;
   invshrink=1./shrink+Input.epsilon;
   shrink=1./invshrink;
  
  invB0(1+first:end,1+first:end)= diag(shrink);
 
  
 if sum(abs(betamean-oldbetmean)<0.0001)==K && count>20
  check=1;
 end
 
 if count==200
    check=1; 
 end
  
end
%% forecats 
Output.pointf=Xout*betamean;
Output.betamean=betamean;
 Output.betamean_sparse=sparsify(X,betamean);
Output.pointf_sparse=Xout*Output.betamean_sparse;

% plot([y X*betamean])
%plot(betamean)


