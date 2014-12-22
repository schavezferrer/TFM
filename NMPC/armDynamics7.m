
clear all 
clc
close all

mdl_armParam7

Ts = 0.5; % seconds
timeSim = 20; % seconds
samples = floor(timeSim / Ts);

ang = [0 0 0 0 0 0 0];
velAng = [0 0 0 0 0 0 0];
accAng = [0 0 0 0 0 0 0];
baseReaction = [0 0 0 0 0 0 0];
jointTorques = [0 0 0 0 0 0 0];

u = zeros(samples,7); %torque controled
% u(:,1) = 0.25;

myData = zeros(samples,21+6+1);

myData(1,1:7) = ang;
myData(1,8:14) = velAng;
myData(1,15:21) = accAng;

for k = 1 : samples
%    velAng = u(k,:);
   
    if(k > 0 && k < 3)
        u(k,2) = 0.05;
      
    elseif(k > 15 && k < 18)
        u(k,2) = -0.05;
    end
        
   c = arm.coriolis(ang, velAng);
%    c*velAng'
   
    accAng = arm.accel(ang,velAng,arm.gravload(ang)+u(k,:)*1 + (c*velAng')');
    accAng(1) = 0;
%    accAng = arm.accel(ang,velAng,arm.gravload(ang));

%    [jointTorques baseReaction] = arm.rne(ang,velAng,accAng');
   
   velAng = velAng + (accAng)'*Ts;
   ang = ang + velAng*Ts;

   [jointTorques baseReaction] = arm.rne(ang,velAng,accAng');

   
   myData(k+1,1:7) = ang;
   myData(k+1,8:14) = velAng;
   myData(k+1,15:21) = accAng';
   myData(k+1,22:27) = baseReaction;
   
   
    myData(k+1,28) = k*Ts;
end

figure(1) 

subplot(3,2,1)
grid on
ylabel('Force Reaction X')
hold on
plot(myData(:,28),myData(:,22),'r')

subplot(3,2,3)
grid on
ylabel('Force Reaction Y')
hold on
plot(myData(:,28),myData(:,23),'g')

subplot(3,2,5)
grid on
ylabel('Force Reaction Z')
hold on
plot(myData(:,28),myData(:,24),'b')

subplot(3,2,2)
grid on
ylabel('Torque Reaction X')
hold on
plot(myData(:,28),myData(:,25),'r')

subplot(3,2,4)
grid on
ylabel('Torque Reaction Y')
hold on
plot(myData(:,28),myData(:,26),'g')

subplot(3,2,6)
grid on
ylabel('Torque Reaction Z')
hold on
plot(myData(:,28),myData(:,27),'b')


figure(2) 

for k = 1 : 7

    for j = 1:3
        
        subplot(7,3,j+(k-1)*3)
        grid on
        hold on
        if j == 1 
           
            plot(myData(:,28),myData(:,k),'r') 
            
        end
        if j == 2
            
            plot(myData(:,28), myData(:,k+7),'g')
            
        end
        if j == 3
            
            plot(myData(:,28),myData(:,k+14),'b')
        
        end

    end
    
end

% 
% 
% hold on
% arm.plot(ang);





