function [ ang, vel,accAng,baseReaction] = armPerturbation( ang, vel, controlSignal, Ts, arm)
%ARMPERTURBATION Summary of this function goes here
%   Detailed explanation goes here

   
   c = arm.coriolis(ang, vel);

   accAng = arm.accel(ang,vel,arm.gravload(ang)+controlSignal+(c*vel')');
%    vel = controlSignal;
%    accAng = arm.accel(ang,vel,arm.gravload(ang));

   vel = vel + (accAng)'*Ts;
   ang = ang + vel*Ts;
   
   [~ , baseReaction] = arm.rne(ang,vel,accAng');


end

