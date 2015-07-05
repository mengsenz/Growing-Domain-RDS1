function [a,m] =assemble(nu,ne,W,Qx,element,node)


m=zeros(nu,nu);       % m matrix
a=sparse(nu,nu);      % a matrix

for e=1:ne
   
    for q = 1:size(W,1)
        w = W(q);
        xpt = Qx(e,q);
        [basis,derivbasis]=lbasis(xpt,e,element,node); 
        
       for iln=1:2
           in = element(e,iln);
        
           
       
               for jln=1:2
                   jn=element(e,jln);
                
                   aji = derivbasis(iln)*derivbasis(jln);
                   a(jn,in)= a(jn,in)+ aji*w;
                   
                   
                   mji = basis(iln)*basis(jln);
                   m(jn,in)= m(jn,in)+ mji*w;

              
               end
       end
    end
end
