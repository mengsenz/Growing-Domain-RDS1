function guv =assemble_guv(nu,ne,W,Qx,element,node,bb,uk,vk,gamma)

guv=zeros(nu,1);       % m matrix
   
for e=1:ne
   
    for q = 1:size(W,1)
        w = W(q);
        xpt = Qx(e,q);
        [basis,~]=lbasis(xpt,e,element,node); 

       for iln=1:2
           in = element(e,iln);
        

                   mji = gamma*(bb - (uk(in)^2)*vk(in));
                   guv(in)= guv(in)+ mji*basis(iln)*w;

         
       end
    end
end
