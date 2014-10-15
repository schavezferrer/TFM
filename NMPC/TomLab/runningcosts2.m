function cost = runningcosts2(u,Prob)
   
%     cost =   (x(3) + 4)^2 + (x(1) - 0.5)^2;
    system2(0, x, u, 0.01)


  cost = (x(1)-1)^2+(x(2)-1)^2;
   
   
