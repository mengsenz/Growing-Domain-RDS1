function [basis,derivbasis]=lbasis(x,e,element,node)


x1 = node(element(e,1));
x2 = node(element(e,2));

d = (x1-x2);


t1 = (x-x2)/d;
t2 = (x1-x)/d;

basis = [t1;
         t2];
       
            
        
 derivbasis = [ 1/d;                  
               -1/d];                
                