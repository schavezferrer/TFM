%clear all 
%clc

mdl_quadcopter

D(:,1) = [quad.d;0;quad.h];          %Di         Rotor hub displacements             1x3
D(:,2) = [0;quad.d;quad.h];
D(:,3) = [-quad.d;0;quad.h];
D(:,4) = [0;-quad.d;quad.h];

pos = [0 0 -4];
ang = [0 0 0];
linVel = [0 0 0];
angVel = [0 0 0];

state = [pos ang linVel angVel]';

Ts = 0.05; % seconds
timeSim = 10; % seconds
samples = floor(timeSim / Ts);

myData = zeros(12,samples);
derivatives = zeros(12,samples);

myData(:,1) = state'; 

%w = ones(samples,4)*300;
w = result(:,14:17);

thrust = zeros(samples,3,4);
Q = zeros(samples,3,4);
tau = zeros(samples,3,4);
a1s = zeros(samples,4);
b1s = zeros(samples,4);

for k = 1 : samples
    
    %UPDATE DYANAMICS
    
    R = [cos(state(5))*cos(state(4))    sin(state(6))*sin(state(5))*cos(state(4))-cos(state(6))*sin(state(4))   cos(state(6))*sin(state(5))*cos(state(4))+sin(state(6))*sin(state(4));   %BBF > Inertial rotation matrix
         cos(state(5))*sin(state(4))    sin(state(6))*sin(state(5))*sin(state(4))+cos(state(6))*cos(state(4))   cos(state(6))*sin(state(5))*sin(state(4))-sin(state(6))*cos(state(4));
         -sin(state(5))                 sin(state(6))*cos(state(5))                                             cos(state(6))*cos(state(5))];

   iW = [0              sin(state(6))                   cos(state(6));             %inverted Wronskian
         0              cos(state(6))*cos(state(5))     -sin(state(6))*cos(state(5));
         cos(state(5))  sin(state(6))*sin(state(5))     cos(state(6))*sin(state(5))] / cos(state(5));
                                         
    for motor = 1:4

        Vr = cross(state(10:12),D(:,motor)) + state(7:9); % Calculate of the rotor (Body Frame)
    
        mu = sqrt(sum(Vr(1:2).^2)) / (abs(w(k,motor))*quad.r);  %Magnitude of mu, planar components
        lc = Vr(3) / (abs(w(k,motor))*quad.r);   %Non-dimensionalised normal inflow
        li = mu; %Non-dimensionalised induced velocity approximation
        j = atan2(Vr(2),Vr(1));  %Sideslip azimuth relative to e1 (zero over nose)
        J = [cos(j) -sin(j);
            sin(j) cos(j)];  %BBF > mu sideslip rotation matrix
        
        %Flapping effect
        
        beta = [((8/3*quad.theta0 + 2*quad.theta1)*mu - 2*(lc)*mu)/(1-mu^2/2); %Longitudinal flapping
            0];%sign(w(i)) * (4/3)*((quad.Ct/quad.sigma)*(2*mu*quad.gamma/3/quad.a)/(1+3*quad.e/2/quad.r) + li)/(1+mu^2/2)]; %Lattitudinal flapping (note sign)
        beta = J'*beta;  %Rotate the beta flapping angles to longitudinal and lateral coordinates.
        o = state(10:12);
        a1s(k,motor) = beta(1) - 16/quad.gamma/abs(w(k,motor)) * o(2);
        b1s(k,motor) = beta(2) - 16/quad.gamma/abs(w(k,motor)) * o(1);

        thrust(k,:,motor) = quad.Ct*quad.rho*quad.A*quad.r^2*w(k,motor)^2* [-cos(b1s(k,motor))*sin(a1s(k,motor)); sin(b1s(k,motor));-cos(a1s(k,motor))*cos(b1s(k,motor))];
       
        Q(k,:,motor) = -quad.Cq*quad.rho*quad.A*quad.r^3*w(k,motor)*abs(w(k,motor)) * [0;0;1];     %Rotor drag torque - note that this preserves w(i) direction sign
        
        tau(k,:,motor) = cross(thrust(k,:,motor),D(:,motor));    %Torque due to rotor thrust
   
    end
    
    %Derivative
    
    dz =  eye(3)*state(7:9);
    
    dn = iW*state(10:12);
    
    totalThrust = thrust(k,:,1)'+ thrust(k,:,2)'+ thrust(k,:,3)'+ thrust(k,:,4)';
    dv = (quad.g*[0;0;1]  + R*(1/quad.M)*totalThrust);
 
    totalTau = tau(k,:,1)'+tau(k,:,2)'+tau(k,:,3)'+tau(k,:,4)';
    totalDrag = Q(k,:,1)'+Q(k,:,2)'+Q(k,:,3)'+Q(k,:,4)';
    do = inv(quad.J)*(cross(-state(10:12),quad.J*state(10:12)) + totalTau + totalDrag);
    
    derivatives(:,k) = [dz' dn' dv' do']';
       
    %UPDATE LINEAR VELOCITIES (Inertial Frame)
    linVel = state(7:9) + dv*Ts; 
    
    %UPDATE ANGLES (Body Frame)
    
    angVel = state(10:12) + do*Ts;
    
    %UPDATE POSITIONS (Inertial Frame)
    pos = state(1:3) + dz*Ts;
    
    %UPDATE ANGLES (Inertial Frame)
    ang = state(4:6) + dn*Ts; 
        
    derivatives(:,k) = [dz' dn' dv' do']';
    
    %POST UPDATE
    state = [pos' ang' linVel' angVel']'; %(I I I B)

    myData(:,k+1) = state;
end

myData = myData';
trueData = result(:,2:13);
derivatives = derivatives';
threshold = 0.0001;

for i = 1:samples

    state = myData(i,:);
    
 
    iW = [0              sin(state(6))                   cos(state(6));
        0               cos(state(6))*cos(state(5))     -sin(state(6))*cos(state(5)); 
        cos(state(5))   sin(state(6))*sin(state(5))     cos(state(6))*sin(state(5))]/cos(state(5));
    
    R = [cos(state(5))*cos(state(4)) sin(state(6))*sin(state(5))*cos(state(4))-cos(state(6))*sin(state(4)) cos(state(6))*sin(state(5))*cos(state(4))+sin(state(6))*sin(state(4));   %BBF > Inertial rotation matrix
         cos(state(5))*sin(state(4)) sin(state(6))*sin(state(5))*sin(state(4))+cos(state(6))*cos(state(4)) cos(state(6))*sin(state(5))*sin(state(4))-sin(state(6))*cos(state(4));
         -sin(state(5))         sin(state(6))*cos(state(5))                            cos(state(6))*cos(state(5))];
   
    myData(i,7:9) =  R'*myData(i,7:9)'; % Lo paso a Body Frame
    myData(i,10:12) = iW*myData(i,10:12)'; % Lo paso a Inertial Frame

end

close all
figure
pMyData = 1;
pTrueData = 1;

subplot(2,2,1)
if pMyData plot3(myData(1:samples,1),myData(1:samples,2),abs(myData(1:samples,3)),'r+'); end
hold on
grid on
if pTrueData plot3(trueData(1:samples,1),trueData(1:samples,2),abs(trueData(1:samples,3)),'b-+'); end
title('Position')
xlabel('Xpos')
ylabel('Ypos')
zlabel('Zpos')

subplot(2,2,2)
if pMyData plot3(myData(1:samples,4),myData(1:samples,5),myData(1:samples,6),'r+'); end
hold on
grid on
if pTrueData plot3(trueData(1:samples,4),trueData(1:samples,5),trueData(1:samples,6),'b-+'); end
title('Angle')
xlabel('Yaw')
ylabel('Pitch')
zlabel('Roll')

subplot(2,2,3)
if pMyData plot3(myData(1:samples,7),myData(1:samples,8),abs(myData(1:samples,9)),'r+'); end
hold on
grid on
if pTrueData plot3(trueData(1:samples,7),trueData(1:samples,8),abs(trueData(1:samples,9)),'b-+'); end
title('linVel')
xlabel('XVel')
ylabel('YVel')
zlabel('ZVel')

subplot(2,2,4)
if pMyData plot3(myData(1:samples,10),myData(1:samples,11),myData(1:samples,12),'r+'); end
hold on
grid on
if pTrueData plot3(trueData(1:samples,10),trueData(1:samples,11),trueData(1:samples,12),'b-+'); end
title('AngVel')
xlabel('YawVel')
ylabel('PitchVel')
zlabel('RollVel')


clear pos ang linVel angVel tout longitudinal transversal J R Tsa1 b1 flapping i k
clear i j l lc li motor mu samples vr a1 timeSim state Vr Ts D Ixx Iyy Izz iW
clear threshold
    