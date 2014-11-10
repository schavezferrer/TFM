function cost = runningcosts(t, x, u, lastU)
   
%     cost =   (x(3) + 4)^2 + (x(1) - 0.5)^2;
    
    cost = (x(1) - 1)^2 + (x(2) - 1)^2 + 1.5*(x(3) + 4)^2; %+ ...
%            (x(4) - 0)^2 + (x(5) - 0)^2 + (x(6) - 0)^2;
   
   
end