function fuv =assemble_fuv(nu,ne,W,Qx,element,node,uk,vk,aa,gamma)

fuv=zeros(nu,1);       % m matrix
  
for e=1:ne
   
    for q = 1:size(W,1)
        w = W(q);
        xpt = Qx(e,q);
        [basis,~]=lbasis(xpt,e,element,node); 

       for iln=1:2
           in = element(e,iln);
        

                   mji = gamma*(aa - uk(in) + (uk(in)^2)*vk(in));
                   fuv(in)= fuv(in)+ mji*basis(iln)*w;

         
       end
    end
end
