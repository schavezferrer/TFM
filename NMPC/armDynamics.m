
clear all 
clc
close all

mdl_armParam

Ts = 0.25; % seconds
timeSim = 10; % seconds
samples = floor(timeSim / Ts);

ang = [0 0 0 0 0 0];
velAng = [0 0 0 0 0 0];
accAng = [0 0 0 0 0 0];
baseReaction = [0 0 0 0 0 0];
jointTorques = [0 0 0 0 0 0];

u = zeros(samples,6); %torque controled
u(1:10,1) = 1;
u(11:20,1) = 0;
u(21:30,1) = 0;
u(31:40,1) = 0;
myData = zeros(samples,18+6+1);

myData(1,1:6) = ang;
myData(1,7:12) = velAng;
myData(1,13:18) = accAng;

for k = 1 : samples

   accAng = arm.accel(ang,velAng,arm.gravload(ang)+u(k,:))
   [jointTorques baseReaction] = arm.rne(ang,velAng,accAng');
   
   velAng = velAng + (accAng)'*Ts;
   ang = ang + velAng*Ts;
   
   myData(k+1,1:6) = ang;
   myData(k+1,7:12) = velAng;
   myData(k+1,13:18) = accAng';
   myData(k+1,19:24) = baseReaction;
   
   
    myData(k+1,25) = k*Ts;
end

figure(1) 

subplot(3,2,1)
grid on
ylabel('Force Reaction X')
hold on
plot(myData(:,25),myData(:,19),'r')

subplot(3,2,3)
grid on
ylabel('Force Reaction Y')
hold on
plot(myData(:,25),myData(:,20),'g')

subplot(3,2,5)
grid on
ylabel('Force Reaction Z')
hold on
plot(myData(:,25),myData(:,21),'b')

subplot(3,2,2)
grid on
ylabel('Torque Reaction X')
hold on
plot(myData(:,25),myData(:,22),'r')

subplot(3,2,4)
grid on
ylabel('Torque Reaction Y')
hold on
plot(myData(:,25),myData(:,23),'g')

subplot(3,2,6)
grid on
ylabel('Torque Reaction Z')
hold on
plot(myData(:,25),myData(:,24),'b')


figure(2) 

for k = 1 : 6

    for j = 1:3
        
        subplot(6,3,j+(k-1)*3)
        grid on
        hold on
        if j == 1 
           
            plot(myData(:,25),myData(:,k),'r') 
            
        end
        if j == 2
            
            plot(myData(:,25), myData(:,k+6),'g')
            
        end
        if j == 3
            
            plot(myData(:,25),myData(:,k+12),'b')
        
        end

    end
    
end



% hold on
% arm.plot(ang);





