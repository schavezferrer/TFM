function [y] = system2(t, x, u, T)
  
    
    mdl_quadcopterParam

    D(:,1) = [quad.d;0;quad.h];          %Di         Rotor hub displacements             1x3
    D(:,2) = [0;quad.d;quad.h];
    D(:,3) = [-quad.d;0;quad.h];
    D(:,4) = [0;-quad.d;quad.h];
   
    [r c] = size(x);
    
    if(c > 1) state = x';
    else state = x;   
    end
   
    w = u;
    
    %UPDATE DYANAMICS
    
    R = [cos(state(5))*cos(state(4))    sin(state(6))*sin(state(5))*cos(state(4))-cos(state(6))*sin(state(4))   cos(state(6))*sin(state(5))*cos(state(4))+sin(state(6))*sin(state(4));   %BBF > Inertial rotation matrix
         cos(state(5))*sin(state(4))    sin(state(6))*sin(state(5))*sin(state(4))+cos(state(6))*cos(state(4))   cos(state(6))*sin(state(5))*sin(state(4))-sin(state(6))*cos(state(4));
         -sin(state(5))                 sin(state(6))*cos(state(5))                                             cos(state(6))*cos(state(5))];

   iW = [0              sin(state(6))                   cos(state(6));             %inverted Wronskian
         0              cos(state(6))*cos(state(5))     -sin(state(6))*cos(state(5));
         cos(state(5))  sin(state(6))*sin(state(5))     cos(state(6))*sin(state(5))] / cos(state(5));
                                         
    for motor = 1:4

        thrust(:,motor) = quad.Ct*quad.rho*quad.A*quad.r^2*w(motor)^2 * [0;0;-1];
        Q(:,motor) = -quad.Cq*quad.rho*quad.A*quad.r^3*w(motor)*abs(w(motor)) * [0;0;1];     %Rotor drag torque - note that this preserves w(i) direction sign
        tau(:,motor) = cross(thrust(:,motor),D(:,motor));    %Torque due to rotor thrust
   
    end
    
    %Derivative
    
    dz =  eye(3)*state(7:9);
    
    dn = iW*state(10:12); % Inertial frame
    
    totalThrust = thrust(:,1) + thrust(:,2) + thrust(:,3) + thrust(:,4);
    dv = (quad.g*[0;0;1]  + R*(1/quad.M)*totalThrust); % Inertial Frame
 
    totalTau = tau(:,1)+tau(:,2)+tau(:,3)+tau(:,4);
    totalDrag = Q(:,1)+Q(:,2)+Q(:,3)+Q(:,4);
    do = inv(quad.J)*(cross(-state(10:12),quad.J*state(10:12)) + totalTau + totalDrag); % Body Frame ?
    
    %UPDATE LINEAR VELOCITIES (Inertial Frame)
    linVel = state(7:9) + dv*T; 
    
    %UPDATE ANGLES (Body Frame)
 
    angVel = state(10:12) + do*T;
    
    %UPDATE POSITIONS (Inertial Frame)
    pos = state(1:3) + dz*T;
    
    %UPDATE ANGLES (Inertial Frame)
    ang = state(4:6) + dn*T; 
        
    %POST UPDATE
    y = [pos' ang' linVel' angVel']; 
    
   

end