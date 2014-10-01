

function [allData, t, x, u] = nmpcSystemSimple

    addpath('./nmpcroutine');
    clear all;
    clc;
    close all;

    mpciterations = 1; % Horizonte de prediccion
    N             = 8; % Horizonte de control
    T             = 0.25; % Tiempo de muestreo
    tmeasure      = 0.0; 
    posIni = [0 0 -4]; % Posición inicial del quadrotor
    angIni = [0 0 0]; % Orientación inicial del quadrotor
    linVelIni = [0 0 0]; % Velocidad lineal inicial del quadrotor
    angVelIni = [0 0 0]; % Velocidad angular inicial del quadrotor
  
    y0 = [0 0 0 0 0 0];
    y1 = [0 0 0 0 0 0];
    
    yAnt = [];
    
    xmeasure = [posIni angIni linVelIni angVelIni];
    k = 855; % Velocidad angular de los motores inicial
    w1 = k*ones(1,N);
    w2 = -k*ones(1,N);
    w3 = k*ones(1,N);
    w4 = -k*ones(1,N);
    u0 = [w1;w2;w3;w4];
    
    iprint        = 5;
    tol_opt       = 1e-8;
    opt_option    = 0;
    type          = 'difference equation';
    atol_ode_real = [];
    rtol_ode_real = [];
    atol_ode_sim  = [];
    rtol_ode_sim  = [];
    
    currTime = 0;
    timeSim = 10;
    currSample = 1;
    totalSamples = floor(timeSim/T);

    allData =  repmat(struct('timePredicted',{},'xPredicted',{}, ...
        'uPredicted',{},'xReal',{},'sampleTime',{},'currTime',{},...
        'fval',{}),totalSamples, 1 );
    lastU = u0(1:4,1);
    while currTime <= timeSim
    % fprintf('--------------------------------------------------\n'); 
    % Aplico el nmpc
    
     [t, x, u, fval] = nmpc2(@runningcosts, @terminalcosts, @constraints, ...
         @terminalconstraints, @linearconstraints, @system, ...
         mpciterations, N, T, tmeasure, xmeasure, u0,lastU, ...
         tol_opt, opt_option, ...
         type, atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim, ...
         iprint, @printHeader, @printClosedloopData, @plotTrajectories)

    % Evoluciono el sistema
   % fprintf('--------------------------------------------------\n'); 
    %fprintf('Evolucionando el Sistema\n');
    y = system(t,x(1,:),u(1:4,1),T);
    
    
    % Actualizo las variables para el nmpc
    %fprintf('--------------------------------------------------\n'); 
    %fprintf('Actualizando variables del nmpc\n');
    allData(currSample).timePredicted = t';
    allData(currSample).xPredicted = x;
    allData(currSample).uPredicted = reshape(u',4,mpciterations)';
    allData(currSample).xReal = y;
    allData(currSample).sampleTime = T;
    allData(currSample).currTime = tmeasure;
    allData(currSample).fval = fval;
    
    tmeasure = currSample*T; %revisar
    xmeasure = y
    
    % Avanzar un paso de control para ayudar al buscador
    
    uu = reshape(u',4,mpciterations)';
    
    w1 = uu(1,1)*ones(1,N);
    w2 = uu(1,2)*ones(1,N);
    w3 = uu(1,3)*ones(1,N);
    w4 = uu(1,4)*ones(1,N);
    
%     w1 = [uu(2:mpciterations,1);uu(mpciterations,1)];
%     w2 = [uu(2:mpciterations,2);uu(mpciterations,2)];
%     w3 = [uu(2:mpciterations,3);uu(mpciterations,3)];
%     w4 = [uu(2:mpciterations,4);uu(mpciterations,4)];
    u0 = [w1;w2;w3;w4];   
    
    
%     k = 800; % Velocidad angular de los motores inicial
%     w1 = k*ones(1,N);
%     w2 = -k*ones(1,N);
%     w3 = k*ones(1,N);
%     w4 = -k*ones(1,N);
%     u0 = [w1;w2;w3;w4];
    
    
    currSample = currSample + 1;
    currTime = tmeasure;
    lastU = u(1:4,1);
    printInfo(y,allData,currSample-1); 
    
    fprintf('--------------------------------------------------\n'); 
    fprintf('Porcentaje Completado: %i\n\n\n',(currSample-2.0)/totalSamples*100);
   end
    
    rmpath('./nmpcroutine');
end




% Function

function cost = runningcosts(t, x, u, lastU)
   
    cost =   (x(3) + 4)^2 + (x(7) - 1.0)^2;
    
end

function cost = terminalcosts(t, x)
   cost = 0;   
end

function [c,ceq] = constraints(t, x, u,N)
    
    k = 100;
    v = 100;
   
    ub = [ 100 100 100 100 100 100 v v 10 k k k];
       
    lb = [-100 -100 -100 -k -k -k -v -v -10 -k -k -k];
    
    [nothing numVar] = size(x);
    
    c = zeros(numVar*2,1);
    x(3) = -1*x(3);
    x(9) = -1*x(9);
  %     for i = 1:2
%         for j = 1:12*2
%             
%             if mod(j,2) c((i-1)*12*2+j) = x((j+1)/2)-ub((j+1)/2);
%             else c((i-1)*12*2+j) = -x(j/2) + lb(j/2);
%             end
%             
%         end
%     end
%     
    
    i = 1;
    for j = 1:numVar*2

        if mod(j,2) c((i-1)*12*2+j) = x((j+1)/2)-ub((j+1)/2);
        else c((i-1)*12*2+j) = -x(j/2) + lb(j/2);
        end

    end
    
    ceq = [];
    
end

function [c,ceq] = terminalconstraints(t, x)
    c   = [];
    ceq = [];
end

function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
    A = [];
    b = []; 
    Aeq = [];
    beq = []; 
    ub = [1000  -200     1000   -200];
    lb = [200    -1000  200     -1000];




end

function [y] = system(t, x, u, T)
  
    
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of output format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printHeader()
   % fprintf('   k  |      u(k)        x(1)        x(2)     Time\n');
    %fprintf('--------------------------------------------------\n');
end

function printClosedloopData(mpciter, u, x, t_Elapsed)
    %fprintf(' %3d  | %+11.6f %+11.6f  %+6.3f', mpciter, u, x, t_Elapsed);
end

function plotTrajectories(dynamic, system, T, t0, x0, u, ...
                          atol_ode, rtol_ode, type)
%     [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, ...
%                                           x0, u, atol_ode, rtol_ode, type);
%     
%     samples =  size(t_intermediate');
%     myData = x_intermediate;
%     
%     for i = 1:samples
% 
%         state = x_intermediate(i,:);
% 
% 
%         iW = [0              sin(state(6))                   cos(state(6));
%             0               cos(state(6))*cos(state(5))     -sin(state(6))*cos(state(5)); 
%             cos(state(5))   sin(state(6))*sin(state(5))     cos(state(6))*sin(state(5))]/cos(state(5));
% 
%         R = [cos(state(5))*cos(state(4)) sin(state(6))*sin(state(5))*cos(state(4))-cos(state(6))*sin(state(4)) cos(state(6))*sin(state(5))*cos(state(4))+sin(state(6))*sin(state(4));   %BBF > Inertial rotation matrix
%              cos(state(5))*sin(state(4)) sin(state(6))*sin(state(5))*sin(state(4))+cos(state(6))*cos(state(4)) cos(state(6))*sin(state(5))*sin(state(4))-sin(state(6))*cos(state(4));
%              -sin(state(5))         sin(state(6))*cos(state(5))                            cos(state(6))*cos(state(5))];
% 
%         myData(i,7:9) =  R'*myData(i,7:9)'; % Lo paso a Body Frame
%         myData(i,10:12) = iW*myData(i,10:12)'; % Lo paso a Inertial Frame
% 
%     end                                 
                                      
                                      
                                      
                                      
%                                       figure(1);
%         title('Closed loop trajectory');
%         xlabel('t');
%         ylabel('x');
%         grid on;
%         hold on;
% %         t_intermediate(:)
% %         size(x_intermediate(:))
%         x1 = zeros(size(t_intermediate'),1);
%         
%       
%         x1(:,1) = x_intermediate(:,1);
%         plot(t_intermediate(:),x1(:),'-ok');
%         axis([0 1 0 2]);
%         axis square;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function printInfo(state,allData,currSample)
    iW = [0              sin(state(6))                   cos(state(6));
           0               cos(state(6))*cos(state(5))     -sin(state(6))*cos(state(5)); 
           cos(state(5))   sin(state(6))*sin(state(5))     cos(state(6))*sin(state(5))]/cos(state(5));

        R = [cos(state(5))*cos(state(4)) sin(state(6))*sin(state(5))*cos(state(4))-cos(state(6))*sin(state(4)) cos(state(6))*sin(state(5))*cos(state(4))+sin(state(6))*sin(state(4));   %BBF > Inertial rotation matrix
             cos(state(5))*sin(state(4)) sin(state(6))*sin(state(5))*sin(state(4))+cos(state(6))*cos(state(4)) cos(state(6))*sin(state(5))*sin(state(4))-sin(state(6))*cos(state(4));
             -sin(state(5))         sin(state(6))*cos(state(5))                            cos(state(6))*cos(state(5))];
        
    myData = state;
    myData(7:9) =  R'*myData(7:9)'; % Lo paso a Body Frame
    myData(10:12) = iW*myData(10:12)'; % Lo paso a Inertial Frame
    
    
    scnsize = get(0,'ScreenSize');
    pos = [0 scnsize(4)/2 scnsize(3)/2.5 scnsize(4)/1.5];
    graphStates(myData,pos);
    pos = [scnsize(3)/2.5 scnsize(4) scnsize(3)/3 scnsize(4)/2];
    graphQuadRotor(myData,pos);
    pos = [0 0 scnsize(3)/2.5 scnsize(4)/1.5];
    graphU(allData,currSample,pos);
    pos = [scnsize(3)/2.5 0 scnsize(3)/3 scnsize(4)/2.5];
    graphFval(allData,currSample,pos);

    
end

function graphStates(myData,pos)
    
    fig = figure(1);
    set(fig,'OuterPosition',pos) 
    subplot(2,2,1)
    hold on
    title('Position')
    xlabel('Xpos')
    ylabel('Ypos')
    zlabel('Zpos')
    grid on
    plot3(myData(1),myData(2),abs(myData(3)),'r-+')
    view(3)
    drawnow

    subplot(2,2,2)
    hold on
    title('Angle')
    xlabel('Yaw')
    ylabel('Pitch')
    zlabel('Roll')
    grid on
    plot3(myData(4),myData(5),myData(6),'r-+')
    view(3)
    drawnow
    
    
    subplot(2,2,3)
    hold on
    title('linVel')
    xlabel('XVel')
    ylabel('YVel')
    zlabel('ZVel')
    grid on
    plot3(myData(7),myData(8),abs(myData(9)),'r-+')
    view(3)
    drawnow
    
    
    subplot(2,2,4)
    hold on
    title('angVel')
    xlabel('YawVel')
    ylabel('PitchVel')
    zlabel('RollVel')
    grid on
    plot3(myData(10),myData(11),abs(myData(12)),'r-+')
    view(3)
    drawnow
    
end

function graphQuadRotor(myData,pos)
    
  
    fig = figure(2);
    set(fig,'OuterPosition',pos) 
    clf;
    view(3)
    grid on
    
    state = myData;    
    
    mdl_quadcopterParam

    D(:,1) = [quad.d;0;0];          %Di         Rotor hub displacements             1x3
    D(:,2) = [0;quad.d;0];
    D(:,3) = [-quad.d;0;0];
    D(:,4) = [0;-quad.d;0];
    
    R = [cos(state(5))*cos(state(4)) sin(state(6))*sin(state(5))*cos(state(4))-cos(state(6))*sin(state(4)) cos(state(6))*sin(state(5))*cos(state(4))+sin(state(6))*sin(state(4));   %BBF > Inertial rotation matrix
             cos(state(5))*sin(state(4)) sin(state(6))*sin(state(5))*sin(state(4))+cos(state(6))*cos(state(4)) cos(state(6))*sin(state(5))*sin(state(4))-sin(state(6))*cos(state(4));
             -sin(state(5))         sin(state(6))*cos(state(5))                            cos(state(6))*cos(state(5))];
    
    hold on
    plot3(myData(1),myData(2),abs(myData(3)),'*b')
    
    p = R*D(:,1) + state(1:3)';
    plot3(p(1),p(2),abs(p(3)),'ok')
    plot3([myData(1) p(1)],[myData(2) p(2)],[abs(myData(3)) abs(p(3))],'k');
    
    p = R*D(:,2) + state(1:3)';
    plot3(p(1),p(2),abs(p(3)),'ok')
    plot3([myData(1) p(1)],[myData(2) p(2)],[abs(myData(3)) abs(p(3))],'k');

    p = R*D(:,3) + state(1:3)';
    plot3(p(1),p(2),abs(p(3)),'ok')
    plot3([myData(1) p(1)],[myData(2) p(2)],[abs(myData(3)) abs(p(3))],'k');
 
    p = R*D(:,4) + state(1:3)';
    plot3(p(1),p(2),abs(p(3)),'ok')
    plot3([myData(1) p(1)],[myData(2) p(2)],[abs(myData(3)) abs(p(3))],'k');
    
    
    xAxis = R*[1 0 0]'*0.1 + myData(1:3)';
    plot3([myData(1) xAxis(1)],[myData(2) xAxis(2)],[abs(myData(3)) abs(xAxis(3))],'r');
   
    yAxis = R*[0 1 0]'*0.1 + myData(1:3)';
    plot3([myData(1) yAxis(1)],[myData(2) yAxis(2)],[abs(myData(3)) abs(yAxis(3))],'g');
    
    zAxis = R*[0 0 1]'*0.1 + myData(1:3)';
    plot3([myData(1) zAxis(1)],[myData(2) zAxis(2)],abs([myData(3) zAxis(3)]),'b');
    drawnow
    
%     workspace = [get(gca,'XLim') get(gca,'YLim') get(gca,'ZLim')];
%     
%     minScale = min(abs(workspace));
%     
%     workspace = [-minScale minScale -minScale minScale -minScale minScale]
    
    mdl_armParam
    
    arm.base = [ [R*[1 0 0]'*0.1 R*[0 1 0]'*0.1  [0 0 0.1]'] [myData(1) myData(2) abs(myData(3))]' ;
                zeros(1,3) 1];
            
    arm.plot(qn,'noshadow','nobase','perspective','magscale',0.05);
    drawnow
    
end


function graphU(allData, currSample,pos)

    uOpt = zeros(currSample,4);
    for i = 1:currSample
        uOpt(i,:) = allData(i).uPredicted(1,:);
    end
    
   
    fig = figure(3);
    set(fig,'OuterPosition',pos) 

    
    subplot(4,1,1)
%     plot(1:currSample, uOpt(:,1),'ro-');
    bar(uOpt(:,1),'b')
    grid on

    subplot(4,1,2)
%     plot(1:currSample, uOpt(:,2),'go-');
    bar(uOpt(:,2),'r')
    grid on
    
    subplot(4,1,3)
%     plot(1:currSample, uOpt(:,3),'ro-');
    bar(uOpt(:,3),'b')
    grid on
    
    subplot(4,1,4)
%     plot(1:currSample, uOpt(:,4),'go-');
    bar(uOpt(:,4),'r')
    grid on
    drawnow
end

function graphFval(allData, currSample,pos)

    fval = zeros(currSample,1);
    for i = 1:currSample
      fval(i) = allData(i).fval;
    end
    
   
    fig = figure(4);
    set(fig,'OuterPosition',pos) 
  
    bar(fval(:,1),'b')
    grid on
    drawnow

end


