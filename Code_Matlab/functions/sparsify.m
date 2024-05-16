function[C_sparse]=sparsify(X,C)
[p,q]=size(C);
C_sparse=C ;
pen=2;
for i =1:q 
    for k =1:p        
        mu=1/(abs(C(k,i))^pen);
        norm=sum(X(:,k).^2);
        if abs(C(k,i))*norm<mu 
            C_sparse(k,i)=0;
        else
            C_sparse(k,i)=sign(C(k,i))*(1/norm)*(abs(C(k,i))*norm-mu) ;
        end       
    end
end
    
    
    
   