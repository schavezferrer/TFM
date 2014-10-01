function [ output_args ] = Untitled(SiSmooth  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    n = 500;
   
    u  =[0 0 0 0];

    yy = [0 0 0 0 0 0];
    yyy = [0 0 0 0 0 0];
    
    x = zeros(n,6);
    
    uu = [0 0 0 0];
    uuu = [0 0 0 0];
    
    yant = [yy(1) yyy(1) yy(2) yyy(2) yy(3) yyy(3) yy(4) yyy(4) yy(5) yyy(5) yy(6) yyy(6)];
    uant = [uu(1) uuu(1) uu(2) uuu(2) uu(3) uuu(3) uu(4) uuu(4)]; 
    
    y = zeros(n,6);
    
    y(1,:) = sim(SiSmooth,u,'init',[yant uant]');
    
    Ts = 0.025;
    
    for i = 3:n
     
        yyy = y(i-2,:);
        yy = y(i-1,:);
        
        uuu = uu;
        uu = u;
        
        yant = [yy(1) yyy(1) yy(2) yyy(2) yy(3) yyy(3) yy(4) yyy(4) yy(5) yyy(5) yy(6) yyy(6)];
        uant = [uu(1) uuu(1) uu(2) uuu(2) uu(3) uuu(3) uu(4) uuu(4)]; 
        
        u = [(0-y(i-1,1))*0.9 (0-y(i-1,2)) 0 0];
        
        y(i,:) = sim(SiSmooth,u,'init',[yant uant]');
        
        x(i,:) = x(i-1,:) + y(i,:)*Ts;
        
    end
    
    
    figure(1)
    
    subplot(3,2,1)
    plot(y(:,1),'r');
        
    subplot(3,2,3)
    plot(y(:,2),'g');

    subplot(3,2,5)
    plot(y(:,3),'b');

%     subplot(3,2,2)
%     plot(y(:,4),'r');
%         
%     subplot(3,2,4)
%     plot(y(:,5),'g');
% 
%     subplot(3,2,6)
%     plot(y(:,6),'b');
    
    
    figure(2)
    
    subplot(3,2,1)
    plot(x(:,1),'r');
        
    subplot(3,2,3)
    plot(x(:,2),'g');

    subplot(3,2,5)
    plot(x(:,3),'b');

%     subplot(3,2,2)
%     plot(x(:,4),'r');
%         
%     subplot(3,2,4)
%     plot(x(:,5),'g');
% 
%     subplot(3,2,6)
%     plot(x(:,6),'b');
    
    
end

