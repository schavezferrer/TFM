function [ ang, vel,accAng,baseReaction] = armPerturbation( ang, vel, controlSignal, Ts, arm)
%ARMPERTURBATION Summary of this function goes here
%   Detailed explanation goes here

   accAng = arm.accel(ang,vel,arm.gravload(ang)+controlSignal);
%    vel = controlSignal;
%    accAng = arm.accel(ang,vel,arm.gravload(ang));

   [~ , baseReaction] = arm.rne(ang,vel,accAng');

   
   vel = vel + (accAng)'*Ts;
   ang = ang + vel*Ts;
   
   

end

