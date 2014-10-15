function [ A, b_L,b_U ] = linearConstraints( lastU )
%LINEARCONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here
   
    maxRate = 100;
    
    A = [1 0 0 0;
        0 1 0 0;
        0 0 1 0;
        0 0 0 1];
    
    b_L = [-maxRate+lastU(1,1) -maxRate+lastU(2,1) -maxRate+lastU(3,1) -maxRate+lastU(4,1)];
    
    
    b_U = [maxRate+lastU(1,1) maxRate+lastU(2,1) maxRate+lastU(3,1) maxRate+lastU(4,1)];
   
    k = 100;
    b = [k+lastU(1,1) k-lastU(1,1) k+lastU(2,1) k-lastU(2,1) k+lastU(3,1) k-lastU(3,1) k+lastU(4,1) k-lastU(4,1)]; 
    
%     b = [];
%     A = [];

    ub = [1000  -200     1000   -200];
    lb = [200    -1000  200     -1000];




end

