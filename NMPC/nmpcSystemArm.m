

function [allData, t, x, u] = nmpcSystemArm

    addpath('./nmpcroutine');
    clear all;
    clc;
    close all;

    mpciterations = 1; % Horizonte de prediccion
    N             = 10; % Horizonte de control
    T             = 0.2; % Tiempo de muestreo
    tmeasure      = 0.0; 
    posIni = [0 0 -4]; % Posición inicial del quadrotor
    angIni = [0 0 0]; % Orientación inicial del quadrotor
    linVelIni = [0 0 0]; % Velocidad lineal inicial del quadrotor
    angVelIni = [0 0 0]; % Velocidad angular inicial del quadrotor
  
    y0 = [0 0 0 0 0 0];
    y1 = [0 0 0 0 0 0];
    
    yAnt = [];
    
    xmeasure = [posIni angIni linVelIni angVelIni];
    k = 900; % Velocidad angular de los motores inicial
%     w1 = 940*ones(1,N);
%     w2 = -870*ones(1,N);
%     w3 = 865*ones(1,N);
%     w4 = -873*ones(1,N);
    
    w1 = k*ones(1,N);
    w2 = -k*ones(1,N);
    w3 = k*ones(1,N);
    w4 = -k*ones(1,N);
    u0 = [w1;w2;w3;w4];
    
    iprint        = 5;
    tol_opt       = 1e-8;
    opt_option    = 1;
    type          = 'difference equation';
    atol_ode_real = [];
    rtol_ode_real = [];
    atol_ode_sim  = [];
    rtol_ode_sim  = [];
    
    currTime = 0;
    timeSim = 10;
    currSample = 1;
    totalSamples = floor(timeSim/T)+N;

    allData =  repmat(struct('timePredicted',{},'xPredicted',{}, ...
        'uPredicted',{},'xReal',{},'sampleTime',{},'currTime',{},...
        'fval',{} , 'armPerturbance', {},'totalTau',{},...
        'armStates',{},'ISE',{}),totalSamples, 1 );
    lastU = u0(1:4,1);
    

    mdl_armParam7;
    
    ang = [0 0 0 0 0 0 0];
    velAng = [0 0 0 0 0 0 0];
    baseReaction = zeros(6,totalSamples+N);
    controlSignal = zeros(7,totalSamples+N);
    
    e1 = zeros(1,7);
    e2 = zeros(1,7);
    e3 = zeros(1,7);
    
    [ang, velAng,accAng,reaction] = armPerturbation( ang, velAng, controlSignal(:,1)', T, arm);

    baseReaction(:,1) = reaction;
    armState = [ang velAng accAng'];

    allData(1).armPerturbance = [baseReaction(1:3,1)/1 baseReaction(4:6,1)];
    allData(1).armStates = armState;
    allTarget = zeros(totalSamples+N,7);
    l = 0;
    target = [0 0 0 0 0 0 0];
   for k = 2 : totalSamples+N+1
        
        if l> 25 && l < 200 
            target = [0 1 1 1 1 1 1]*0.0;
            allTarget(k,:) = target;
        elseif l > 200
            target = [0 2 2 2 2 2 2]*0;
        end
        
%         if l> 30 
%             l = 0;
%         end
        
%         target = [0 1 1 1 1 1 1]*0.1*sin(5*k*2*3.1415/totalSamples);
%         target = [0 1 1 1 1 1 1]*2*k/totalSamples;

        e3 = (target - ang) - e1;
        e1 = (target - ang);
        e2 = e2 + e1*T;
            
            
%             vel = 20*e1 + 0*e2 + 40*e3/T;
            
        vel =   0.1*[0 30 30 30 20 20 70].*e1 ... 
              + 0.1*[0 0 0 0 0 0 0].*e2 ...
              + 0.1*[0 50 90 40 15 15 5].*e3/T;


        newSignal =  0.2*[0 0.23 0.24 0.12 0.0038 0.00033 0.00022].*(vel - velAng)/T;


        maxVal = [0 1 1 1 1 1 1]*10;

        for i = 2:7

            if ( newSignal(1,i) > maxVal(1,i)) newSignal(1,i) = maxVal(1,i);
            elseif ( newSignal(1,i) < -maxVal(1,i)) newSignal(1,i) = -maxVal(1,i);
            end
            controlSignal(i,k) = newSignal(1,i);

        end
        l = l+1;

        [ang, velAng,accAng,reaction] = armPerturbation( ang, velAng, controlSignal(:,k)', T, arm);
   
        baseReaction(:,k) = reaction;
        armState = [ang velAng accAng'];
        
        allData(k).armPerturbance = [baseReaction(1:3,k)/1 baseReaction(4:6,k)];
        allData(k).armStates = armState;
        
               
   end
   
   scnsize = get(0,'ScreenSize');
   
   pos = [scnsize(3)*0.5 scnsize(4) scnsize(3)/2.5 scnsize(4)/1.5];
     
   graphArmReaction(allData,totalSamples,pos,6,N,T);
   
   graphArmStates(allData,totalSamples,pos,7,allTarget,N,T);
    t = 0;
    x = xmeasure;
    u = u0';
    y = xmeasure;
    fval = runningcosts(t, x, u, u)
    while currTime <= timeSim
    % fprintf('--------------------------------------------------\n'); 
    % Aplico el nmpc
    
%     uArm = zeros(1,6); %torque controled
%     if(currTime <= 1) uArm(1,1) = 0.1; end 
%     [ang, velAng, baseReaction] = armPerturbation(ang,velAng,uArm,T,arm);
    allData(currSample).timePredicted = t';
    allData(currSample).xPredicted = x;
    if currSample > 1 
        allData(currSample).uPredicted = reshape(u',4,mpciterations)';
    else
         allData(currSample).uPredicted = u;
    end
    
    allData(currSample).xReal = y;
    allData(currSample).sampleTime = T;
    allData(currSample).currTime = tmeasure;
    
%     if currSample > 1 &&  abs(allData(currSample-1).fval - fval)/allData(currSample-1).fval > 200 
%          if allData(currSample-1).fval - fval < 0
%           allData(currSample).fval =   1.2*allData(currSample-1).fval;
%          else
%            allData(currSample).fval =  fval;
%          end
%      
%     else 
        
            allData(currSample).fval = fval;

%     end
    
%    allData(currSample).totalTau = totalTau;
    if currSample > 1 
        allData(currSample).ISE = allData(currSample).fval + allData(currSample-1).ISE;
    else
        allData(currSample).ISE = allData(currSample).fval;
    end
    
     [t, x, u, fval] = nmpcArm(@runningcosts, @terminalcosts, @constraints, ...
         @terminalconstraints, @linearconstraints, @realSystem, ...
         mpciterations, N, T, tmeasure, xmeasure, u0,lastU,baseReaction, ...
         tol_opt, opt_option, ...
         type, atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim, ...
         iprint, @printHeader, @printClosedloopData, @plotTrajectories);

    % Evoluciono el sistema
   % fprintf('--------------------------------------------------\n'); 
    %fprintf('Evolucionando el Sistema\n');
    [y,totalTau] = mdlSystem(t,x(1,:),u(1:4,1),T,baseReaction);
    
    
    % Actualizo las variables para el nmpc
    %fprintf('--------------------------------------------------\n'); 
    %fprintf('Actualizando variables del nmpc\n');
%     allData(currSample).timePredicted = t';
%     allData(currSample).xPredicted = x;
%     allData(currSample).uPredicted = reshape(u',4,mpciterations)';
%     allData(currSample).xReal = y;
%     allData(currSample).sampleTime = T;
%     allData(currSample).currTime = tmeasure;
%     
%     if currSample > 1 &&  abs(allData(currSample-1).fval - fval)/allData(currSample-1).fval > 8 
%          if allData(currSample-1).fval - fval < 0
%           allData(currSample).fval =   1.2*allData(currSample-1).fval
%          else
%            allData(currSample).fval =  fval;
%          end
%      
%     else 
%             allData(currSample).fval = fval;
% 
%     end
%     
%     allData(currSample).totalTau = totalTau;
%     if currSample > 1 
%         allData(currSample).ISE = allData(currSample).fval + allData(currSample-1).ISE;
%     else
%         allData(currSample).ISE = allData(currSample).fval;
%     end

    tmeasure = currSample*T; %revisar
    xmeasure = y;
    
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
    
    if currSample + N-1  == totalSamples
        currTime = timeSim+1;
    end
    
    printInfo(y,allData,currSample-1,T); 
%     
%     uArm = zeros(1,6); %torque controled
%     uArm(1,1) = 0; 
%     
%     [ang, velAng, baseReaction] = armPerturbation(ang,velAng,uArm,T,arm);
    
%     accAng = arm.accel(ang,velAng,arm.gravload(ang)+uArm(1,:));
%     [jointTorques baseReaction] = arm.rne(ang,velAng,accAng');

    
    
    
    
    fprintf('--------------------------------------------------\n'); 
    fprintf('Porcentaje Completado: %i\n\n\n',(currSample-2.0)/totalSamples*100);
   end
    
    rmpath('./nmpcroutine');
end




function cost = terminalcosts(t, x)
   cost = 0;   
end

function [c,ceq] = constraints(t, x, u,N)
    
    k = 0.1;
    v = 1;
   
    ub = [ 5 5 5 k*10 k k v v 5 1 1 1 ];
       
    lb = [-5 -5 -5 -k*10 -k -k -v -v -5 -1 -1 -1];
    
    [nothing numVar] = size(x);
    
    c = zeros(numVar*2,1);
    x(3) = -1*x(3);
    x(9) = -1*x(9);
  
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

function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u, lastU)
  

    A = [1 0 0 0;
         -1 0 0 0
         0 1 0 0;
         0 -1 0 0;
         0 0 1 0;
         0 0 -1 0;
         0 0 0 1;
         0 0 0 -1];
   
    k = 300;
    b = [k+lastU(1,1) k-lastU(1,1) k+lastU(2,1) k-lastU(2,1) k+lastU(3,1) k-lastU(3,1) k+lastU(4,1) k-lastU(4,1)]; 
    
%     b = [];
%     A = [];
    min = 200;
    max = 2000;
    Aeq = [];
    beq = []; 
    ub = [max  -min  max  -min];
    lb = [min  -max  min  -max];




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

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function printInfo(state,allData,currSample,Ts)
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
    
    pos = [0 0 scnsize(3)/2.5 scnsize(4)/1.5];
    graphU(allData,currSample,pos,1,Ts);
    
    pos = [scnsize(3)*0.5 0 scnsize(3)/3 scnsize(4)/2.5];
    graphFval(allData,currSample,pos,2,Ts);
    
    pos = [scnsize(3)*0 scnsize(4) scnsize(3)/2.5 scnsize(4)/1.5];
    graphIndividualStates(allData,currSample,pos,3,Ts);
    
%     pos = [scnsize(3)*0.5 scnsize(4) scnsize(3)/2.5 scnsize(4)/1.5];
%     graphArmReaction(allData,currSample,pos,4);
    
    pos = [scnsize(3)*0.5 0 scnsize(3)/3 scnsize(4)/2.5];
    %graphTotalTau(allData,currSample,pos,5);
    
    graphISE(allData,currSample,8,Ts);
%     graphQuadRotor(myData,pos,6)
%     graphStates(myData,pos,4);
end

function graphIndividualStates(allData,currSample,pos,numFig,Ts)
    
    myData = zeros(currSample,12);
    for i = 1:currSample
       myData(i,:) = allData(i).xReal(1,:);
    end

    fig = figure(numFig);
%     set(fig,'OuterPosition',pos) 
    
    subplot(6,2,1)
    hold on
    xlabel('Time')
    ylabel('Position X')
    grid on
    plot([1:currSample]*Ts,myData(:,1),'r+-')
     axis tight
    drawnow
    
    subplot(6,2,3)
    hold on
    xlabel('Time')
    ylabel('Position Y')
    grid on
    plot([1:currSample]*Ts,myData(:,2),'g+-')
     axis tight
    drawnow
    
    subplot(6,2,5)
    hold on
    xlabel('Time')
    ylabel('Position Z')
    grid on
    plot([1:currSample]*Ts,-myData(:,3),'b+-')
     axis tight
    drawnow
    
    subplot(6,2,2)
    hold on
    xlabel('Time')
    ylabel('Linear Velocity X')
    grid on
    plot([1:currSample]*Ts,myData(:,7),'r+-')
     axis tight
    drawnow
    
    subplot(6,2,4)
    hold on
    xlabel('Time')
    ylabel('Linear Velocity Y')
    grid on
    plot([1:currSample]*Ts,myData(:,8),'g+-')
     axis tight
    drawnow
    
    subplot(6,2,6)
    hold on
    xlabel('Time')
    ylabel('Linear Velocity Z')
    grid on
    plot([1:currSample]*Ts,-myData(:,9),'b+-')
     axis tight
    drawnow
    
    subplot(6,2,7)
    hold on
    xlabel('Time')
    ylabel('Pitch')
    grid on
    plot([1:currSample]*Ts,myData(:,5),'r+-')
     axis tight
    drawnow
    
    subplot(6,2,9)
    hold on
    xlabel('Time')
    ylabel('Roll')
    grid on
    plot([1:currSample]*Ts,myData(:,6),'g+-')
     axis tight
    drawnow
    
    subplot(6,2,11)
    hold on
    xlabel('Time')
    ylabel('Yaw')
    grid on
    plot([1:currSample]*Ts,myData(:,4),'b+-')
     axis tight
    drawnow
    
    subplot(6,2,8)
    hold on
    xlabel('Time')
    ylabel('Angular Velocity X')
    grid on
    plot([1:currSample]*Ts,myData(:,10),'r+-')
     axis tight
    drawnow
    
    subplot(6,2,10)
    hold on
    xlabel('Time')
    ylabel('Angular Velocity Y')
    grid on
    plot([1:currSample]*Ts,myData(:,11),'g+-')
     axis tight
    drawnow
    
    subplot(6,2,12)
    hold on
    xlabel('Time')
    ylabel('Angular Velocity Z')
    grid on
    plot([1:currSample]*Ts,myData(:,12),'b+-')
     axis tight
    drawnow
    
    
   
end

function graphStates(myData,pos,numFig)
    
    fig = figure(numFig);
%     set(fig,'OuterPosition',pos) 
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
    xlabel('Omega X')
    ylabel('Omega Y')
    zlabel('Omega Z')
    grid on
    plot3(myData(10),myData(11),abs(myData(12)),'r-+')
    view(3)
    drawnow
    
   
end

function graphQuadRotor(myData,pos,numFig)
    
  
    fig = figure(numFig);
%     set(fig,'OuterPosition',pos) 
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

function graphU(allData, currSample,pos,numFig,Ts)

    uOpt = zeros(currSample,8);
  
    for i = 1:currSample
        uOpt(i,1:4) = allData(i).uPredicted(1,:);
%         uOpt(i,5:8) = allData(i).xReal(1,13:16);
    end
    
   
    fig = figure(numFig);
%     set(fig,'OuterPosition',pos) 

  
    subplot(4,1,1)
%     plot(1:currSample, uOpt(:,1),'ro-');

  
    bar([1:currSample]*Ts,uOpt(:,1),'b')
     axis tight
    xlabel('Time')
    ylabel({'Rotor 1','Angular Rate(rad/s)'})
    grid on

    subplot(4,1,2)
%     plot(1:currSample, uOpt(:,2),'go-');
   

    bar([1:currSample]*Ts,uOpt(:,2),'r')
     axis tight
     xlabel('Time')
    ylabel({'Rotor 2','Angular Rate(rad/s)'})
    grid on
    
    subplot(4,1,3)
     %     plot(1:currSample, uOpt(:,3),'ro-');
    bar([1:currSample]*Ts,uOpt(:,3),'b')
     xlabel('Time')
      axis tight
    ylabel({'Rotor 3','Angular Rate(rad/s)'})
    grid on
    
    subplot(4,1,4)

%     plot(1:currSample, uOpt(:,4),'go-');
    bar([1:currSample]*Ts,uOpt(:,4),'r')
        xlabel('Time')
         axis tight
    ylabel({'Rotor 4','Angular Rate(rad/s)'})
    grid on
    drawnow
    

end

function graphArmReaction(allData, currSample,pos,numFig,N,Ts)

    armPerturbance = zeros(currSample,8);
    for i = 1:currSample
        armPerturbance(i,1:6) = allData(i).armPerturbance(:);
    end
    
   
    fig = figure(numFig);
    
%     set(fig,'OuterPosition',pos) 
    subplot(3,2,1)
    grid on
    
    xlabel('Time')

    ylabel('Force Reaction X')
    hold on
    plot([1:currSample-N]*Ts,armPerturbance(1:currSample-N,1),'r')
    axis tight
    
    
    subplot(3,2,3)
    grid on
        xlabel('Time')

    ylabel('Force Reaction Y')
    hold on
    plot([1:currSample-N]*Ts,armPerturbance(1:currSample-N,2),'g')
     axis tight
     
    subplot(3,2,5)
    grid on
        xlabel('Time')

    ylabel('Force Reaction Z')
    hold on
    plot([1:currSample-N]*Ts,armPerturbance(1:currSample-N,3),'b')
 axis tight
 
    subplot(3,2,2)
    grid on
        xlabel('Time')

    ylabel('Torque Reaction X')
    hold on
    plot([1:currSample-N]*Ts,armPerturbance(1:currSample-N,4),'r')
 axis tight
 
    subplot(3,2,4)
    grid on
    xlabel('Time')
    ylabel('Torque Reaction Y')
    hold on
    plot([1:currSample-N]*Ts,armPerturbance(1:currSample-N,5),'g')
 axis tight
 
    subplot(3,2,6)
    grid on
    xlabel('Time')
    ylabel('Torque Reaction Z')
    hold on
    plot([1:currSample-N]*Ts,armPerturbance(1:currSample-N,6),'b')
     axis tight

end

function graphTotalTau(allData, currSample,pos,numFig)

    totalTau = zeros(currSample,8);
    for i = 1:currSample
        totalTau(i,1:3) = allData(i).totalTau(:);
    end
    
   
    fig = figure(numFig);
    
     
    subplot(3,1,1)
    grid on
    ylabel('Torque Reaction X')
    hold on
    plot(totalTau(:,1),'r')
 axis tight
    subplot(3,1,2)
    grid on
    ylabel('Torque Reaction Y')
    hold on
    plot(totalTau(:,2),'g')
 axis tight
    subplot(3,1,3)
    grid on
    ylabel('Torque Reaction Z')
    hold on
    plot(totalTau(:,3),'b')
 axis tight
  
    drawnow

end

function graphFval(allData, currSample,pos,numFig,Ts)

    fval = zeros(currSample,1);
    for i = 1:currSample
      fval(i) = allData(i).fval;
    end
    
   
    fig = figure(numFig);
%     set(fig,'OuterPosition',pos) 
   
    bar([1:currSample]*Ts,fval(:,1),'b')
    axis tight
    xlabel('Time')
    ylabel('Cost Function Value')
    grid on
    drawnow

end

function graphArmStates(allData, currSample,pos,numFig,target,N,Ts)

    armStates = zeros(currSample,21);
    
    for i = 1:currSample
      armStates(i,:) = allData(i).armStates;
    end
    
     figure(20)

     for k = 2 : 7
        subplot(6,1,k-1)
        plot([1:currSample-N]*Ts,armStates(1:currSample-N,k),'r')
        hold on 
        plot([1:currSample-N]*Ts,target(1:currSample-N,k),'k');
         axis tight
        legend('Real Value','Reference')
         if k == 2
            title('Angle (rad)')
         end
      
        xlabel('Time');
        ylabel({['Joint ' int2str(k-1)],'radians'});
        grid on
  
     end
     
     figure(21)

     for k = 2 : 7
         subplot(6,1,k-1)
        plot([1:currSample-N]*Ts,armStates(1:currSample-N,k+7),'g')
 axis tight
         if k == 2
            title('Angular Velocity (rad/seconds)')
         end
      
        xlabel('Time');
        ylabel({['Joint ' int2str(k-1)],'rad/seconds'});
         grid on

        
     end
     
    figure(22)

     for k = 2 : 7
         subplot(6,1,k-1)
        plot([1:currSample-N]*Ts,armStates(1:currSample-N,k+14),'b')
 axis tight
         if k == 2
            title('Angular Acceleration(rad/seconds^2)')
         end
      
        xlabel('Time');
        ylabel({['Joint ' int2str(k-1)],'radians/seconds^2'});
        grid on
     end
     
    drawnow
end
    
function graphISE(allData, currSample,numFig,Ts)

    ISE = zeros(currSample,1);
    
    for i = 1:currSample
      ISE(i,1) = allData(i).ISE;
    end
    figure(numFig)
    
    bar([1:currSample]*Ts,ISE(:,1))
    axis tight
    title('Integral of the Square of the Error')
    xlabel('Time');
    ylabel('ISE');
    grid on
    drawnow
    
end
