
clear all 
clc
close all

mdl_armParam

Ts = 0.05; % seconds
timeSim = 1; % seconds
samples = floor(timeSim / Ts);

ang = [0 0 0 0 0 0];
velAng = [1 0 0 0 0 0];
accAng = [0 0 0 0 0 0];
baseReaction = [0 0 0 0 0 0];
jointTorques = [0 0 0 0 0 0];

u = [ones(samples,1)*1 zeros(samples,5)]; %torque controled
myData = zeros(samples,18);

myData(1,1:6) = ang;
myData(1,7:12) = velAng;
myData(1,13:18) = accAng;

for k = 1 : samples
   
   accAng = arm.accel(ang,velAng,u(k,:));
   [jointTorques baseReaction] = arm.rne(ang,velAng,accAng');
   
   velAng = velAng + (accAng)'*Ts;
   ang = ang + velAng*Ts;
   
   myData(k+1,1:6) = ang;
   myData(k+1,7:12) = velAng;
   myData(k+1,13:18) = accAng';

  
end

figure(1) 

hold on
arm.plot(ang);
plot(myData(:,1));
