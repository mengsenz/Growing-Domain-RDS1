function [W,Qx] = quadrature3(element,node,ne)

Qx = zeros (ne,3);
W = zeros(3,1);
 h = (node(2)-node(1));
      
      W(1)= 0.2778*h;   %5*h/18;
      W(2)= 0.4444*h;  %8h/18
      W(3)=W(1);
        
        for e=1:ne
          x1 = node(element(e,1),1);
          Qx(e,1)=x1+0.1127*h;      % 0.5-sqrt(0.15)
          Qx(e,2)=x1+0.5*h;         % 0.5+sqrt(0.15)  
          Qx(e,3)=x1+0.8873*h;      
      end
      
      
          
          
          
      
        