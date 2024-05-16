function Xtilde= cholkernel(X)

%X=randn(10,20);

[n,~]=size(X);

Kernel=ones(n,n);

for i=1:n
for k=i:n
Kernel(i,k)=(X(i,:)-X(k,:))*(X(i,:)-X(k,:))';
Kernel(k,i)=Kernel(i,k);
end
end

mask = tril(true(size(Kernel)),-1);
kappa=median(1./Kernel(mask));

kernel=exp(-kappa*0.5*Kernel);

Xtilde=chol(kernel)'; %Xtilde*Xtilde'=kernel->Xtilde ist untere dreiecksmatrix
