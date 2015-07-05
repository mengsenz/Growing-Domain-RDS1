function [fuv, guv] = assemble_reaction(nu,ne,W,Qx,element,node,uk,vk,reactionFun) %, aa,gamma

fuv=zeros(nu,1);       % activator
guv=zeros(nu,1);       % inhibitor
for e=1:ne
   
    for q = 1:size(W,1)
        w = W(q);
        xpt = Qx(e,q);
        [basis,~]=lbasis(xpt,e,element,node); 

       for iln=1:2
           in = element(e,iln);

           mji_fg = reactionFun(uk(in),vk(in)) ; 
           fuv(in)= fuv(in)+ mji_fg(1)*basis(iln)*w;
           guv(in)= guv(in)+ mji_fg(2)*basis(iln)*w;
         
       end
    end
end
