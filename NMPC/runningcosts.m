function cost = runningcosts(t, x, u, lastU)
   
%     cost =   (x(3) + 4)^2 + (x(1) - 0.5)^2;
    
    cost = (x(1) - 0.5)^2 + (x(2) - 0.5)^2 + 1.5*(x(3) + 4)^2;
   
   
end