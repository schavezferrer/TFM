

function [allData, t, x, u] = nmpcSystemArm

    addpath('./nmpcroutine');
    clear all;
    clc;
    close all;

    mpciterations = 1; % Horizonte de prediccion
    N             = 8; % Horizonte de control
    T             = 3/N; % Tiempo de muestreo
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
    opt_option    = 2;
    type          = 'difference equation';
    atol_ode_real = [];
    rtol_ode_real = [];
    atol_ode_sim  = [];
    rtol_ode_sim  = [];
    
    currTime = 0;
    timeSim = 20;
    currSample = 1;
    totalSamples = floor(timeSim/T);

    allData =  repmat(struct('timePredicted',{},'xPredicted',{}, ...
        'uPredicted',{},'xReal',{},'sampleTime',{},'currTime',{},...
        'fval',{} , 'armPerturbance', {},'totalTau',{},...
        'armStates',{}),totalSamples, 1 );
    lastU = u0(1:4,1);
    

    mdl_armParam;
    
    ang = [0 0 0 0 0 0];
    velAng = [0 0 0 0 0 0];
    baseReaction = zeros(6,totalSamples+N);
    controlSignal = zeros(6,totalSamples+N);
    
    e1 = zeros(1,6);
    e2 = zeros(1,6);
    e3 = zeros(1,6);
    
    [ang, velAng,accAng,reaction] = armPerturbation( ang, velAng, controlSignal(:,1)', T, arm);

    baseReaction(:,1) = reaction;
    armState = [ang velAng accAng'];

    allData(1).armPerturbance = [baseReaction(1:3,1)/1 baseReaction(4:6,1)];
    allData(1).armStates = armState;

    
   for k = 2 : totalSamples+N
   
        target = [1 1 1 1 1 1]*1;
        e3 = (target - ang) - e1;
        e1 = (target - ang);
        e2 = e2 + e1*T;
            
            
%             vel = 20*e1 + 0*e2 + 40*e3/T;
            
        vel = [20 15 10 4 1 0.9].*e1 + [40 35 30 1 1 1].*e3/T;


        newSignal = velAng*0 + [0.06 0.06 0.06 0.00015 0.0003 0.0002].*(vel - velAng);


        maxVal = [0.75 0.75 0.75 0.75 0.75 0.75];

        for i = 1:6

            if ( newSignal(1,i) > maxVal(1,i)) newSignal(1,i) = maxVal(1,i);
            elseif ( newSignal(1,i) < -maxVal(1,i)) newSignal(1,i) = -maxVal(1,i);
            end
            controlSignal(i,k) = newSignal(1,i);

        end
%             controlSignal(:,k)
%             i = 2;
%             if ( newSignal(i) > maxVal(i)) newSignal(i) = maxVal(i);
%             elseif ( newSignal(i) < -maxVal(i)) newSignal(i) = -maxVal(i);
%             end
%             controlSignal(i,k) = newSignal(i);
            
                       
            
%         elseif(k >= 15)
%             controlSignal(1,k) = (0 -  ang(1))*0.001 + 5*(0 -  ang(1))*T;
       
        
%         
%         
%         
        [ang, velAng,accAng,reaction] = armPerturbation( ang, velAng, controlSignal(:,k)', T, arm);
   
        baseReaction(:,k) = reaction;
        armState = [ang velAng accAng'];
        
        allData(k).armPerturbance = [baseReaction(1:3,k)/1 baseReaction(4:6,k)];
        allData(k).armStates = armState;
        
        
       
   end
   
   scnsize = get(0,'ScreenSize');
   
   pos = [scnsize(3)*0.5 scnsize(4) scnsize(3)/2.5 scnsize(4)/1.5];
     
   graphArmReaction(allData,totalSamples,pos,6);
   
   graphArmStates(allData,totalSamples,pos,7);
    
   
   
    while currTime <= timeSim
    % fprintf('--------------------------------------------------\n'); 
    % Aplico el nmpc
    
%     uArm = zeros(1,6); %torque controled
%     if(currTime <= 1) uArm(1,1) = 0.1; end 
%     [ang, velAng, baseReaction] = armPerturbation(ang,velAng,uArm,T,arm);

    
     [t, x, u, fval] = nmpcArm(@runningcosts, @terminalcosts, @constraints, ...
         @terminalconstraints, @linearconstraints, @mdlSystem, ...
         mpciterations, N, T, tmeasure, xmeasure, u0,lastU,baseReaction, ...
         tol_opt, opt_option, ...
         type, atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim, ...
         iprint, @printHeader, @printClosedloopData, @plotTrajectories)

    % Evoluciono el sistema
   % fprintf('--------------------------------------------------\n'); 
    %fprintf('Evolucionando el Sistema\n');
    [y,totalTau] = mdlSystem(t,x(1,:),u(1:4,1),T,baseReaction);
    
    
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
    allData(currSample).totalTau = totalTau;

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
    printInfo(y,allData,currSample-1); 
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




% Function

% function cost = runningcosts(t, x, u, lastU)
%    
%     cost =   (x(3) + 4)^2 + (x(1) - 0.5)^2;
%     
%     cost = (x(1) - 1)^2 + (x(2) - 1)^2 + 1.5*(x(3) + 4)^2;
%    
%    
% end

function cost = terminalcosts(t, x)
   cost = 0;   
end

function [c,ceq] = constraints(t, x, u,N)
    
    k = 0.1;
    v = 1;
   
    ub = [ 5 5 5 k*10 k k v v 10 1 1 1 ];
       
    lb = [-5 -5 -5 -k*10 -k -k -v -v -10 -1 -1 -1];
    
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
    min = 600;
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
    
    pos = [0 0 scnsize(3)/2.5 scnsize(4)/1.5];
    graphU(allData,currSample,pos,1);
    
    pos = [scnsize(3)*0.5 0 scnsize(3)/3 scnsize(4)/2.5];
    graphFval(allData,currSample,pos,2);
    
    pos = [scnsize(3)*0 scnsize(4) scnsize(3)/2.5 scnsize(4)/1.5];
    graphIndividualStates(allData,currSample,pos,3);
    
%     pos = [scnsize(3)*0.5 scnsize(4) scnsize(3)/2.5 scnsize(4)/1.5];
%     graphArmReaction(allData,currSample,pos,4);
    
    pos = [scnsize(3)*0.5 0 scnsize(3)/3 scnsize(4)/2.5];
    graphTotalTau(allData,currSample,pos,5);
    
%     graphQuadRotor(myData,pos,6)
%     graphStates(myData,pos,4);
end

function graphIndividualStates(allData,currSample,pos,numFig)
    
    myData = zeros(currSample,12);
    for i = 1:currSample
       myData(i,:) = allData(i).xReal(1,:);
    end

    fig = figure(numFig);
    set(fig,'OuterPosition',pos) 
   
    subplot(6,2,1)
    hold on
    xlabel('Sample')
    ylabel('Xpos')
    grid on
    plot(myData(:,1),'r+-')
    drawnow
    
    subplot(6,2,3)
    hold on
    xlabel('Sample')
    ylabel('Ypos')
    grid on
    plot(myData(:,2),'g+-')
    drawnow
    
    subplot(6,2,5)
    hold on
    xlabel('Sample')
    ylabel('Zpos')
    grid on
    plot(myData(:,3),'b+-')
    drawnow
    
    subplot(6,2,2)
    hold on
    xlabel('Sample')
    ylabel('XVel')
    grid on
    plot(myData(:,7),'r+-')
    drawnow
    
    subplot(6,2,4)
    hold on
    xlabel('Sample')
    ylabel('YVel')
    grid on
    plot(myData(:,8),'g+-')
    drawnow
    
    subplot(6,2,6)
    hold on
    xlabel('Sample')
    ylabel('ZVel')
    grid on
    plot(myData(:,9),'b+-')
    drawnow
    
    subplot(6,2,7)
    hold on
    xlabel('Sample')
    ylabel('Pitch')
    grid on
    plot(myData(:,5),'r+-')
    drawnow
    
    subplot(6,2,9)
    hold on
    xlabel('Sample')
    ylabel('Roll')
    grid on
    plot(myData(:,6),'g+-')
    drawnow
    
    subplot(6,2,11)
    hold on
    xlabel('Sample')
    ylabel('Yaw')
    grid on
    plot(myData(:,4),'b+-')
    drawnow
    
    subplot(6,2,8)
    hold on
    xlabel('Sample')
    ylabel('PitchVel')
    grid on
    plot(myData(:,10),'r+-')
    drawnow
    
    subplot(6,2,10)
    hold on
    xlabel('Sample')
    ylabel('RollVel')
    grid on
    plot(myData(:,11),'g+-')
    drawnow
    
    subplot(6,2,12)
    hold on
    xlabel('Sample')
    ylabel('YawVel')
    grid on
    plot(myData(:,12),'b+-')
    drawnow
    
    
   
end

function graphStates(myData,pos,numFig)
    
    fig = figure(numFig);
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

function graphQuadRotor(myData,pos,numFig)
    
  
    fig = figure(numFig);
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

function graphU(allData, currSample,pos,numFig)

    uOpt = zeros(currSample,8);
    for i = 1:currSample
        uOpt(i,1:4) = allData(i).uPredicted(1,:);
%         uOpt(i,5:8) = allData(i).xReal(1,13:16);
    end
    
   
    fig = figure(numFig);
    set(fig,'OuterPosition',pos) 

    
    subplot(4,2,1)
%     plot(1:currSample, uOpt(:,1),'ro-');
    bar(uOpt(:,1),'b')
    grid on

    subplot(4,2,3)
%     plot(1:currSample, uOpt(:,2),'go-');
    bar(uOpt(:,2),'r')
    grid on
    
    subplot(4,2,5)
%     plot(1:currSample, uOpt(:,3),'ro-');
    bar(uOpt(:,3),'b')
    grid on
    
    subplot(4,2,7)
%     plot(1:currSample, uOpt(:,4),'go-');
    bar(uOpt(:,4),'r')
    grid on
    drawnow
    
%      subplot(4,2,2)
% %     plot(1:currSample, uOpt(:,1),'ro-');
%     bar(uOpt(:,5),'b')
%     grid on

%     subplot(4,2,4)
% %     plot(1:currSample, uOpt(:,2),'go-');
%     bar(uOpt(:,6),'r')
%     grid on
%     
%     subplot(4,2,6)
%     plot(1:currSample, uOpt(:,3),'ro-');
%     bar(uOpt(:,7),'b')
%     grid on
    
%     subplot(4,2,8)
%     plot(1:currSample, uOpt(:,4),'go-');
%     bar(uOpt(:,8),'r')
%     grid on
%     drawnow
end

function graphArmReaction(allData, currSample,pos,numFig)

    armPerturbance = zeros(currSample,8);
    for i = 1:currSample
        armPerturbance(i,1:6) = allData(i).armPerturbance(:);
    end
    
   
    fig = figure(numFig);
    
    set(fig,'OuterPosition',pos) 
    
    subplot(3,2,1)
    grid on
    ylabel('Force Reaction X')
    hold on
    plot(armPerturbance(:,1),'r')

    subplot(3,2,3)
    grid on
    ylabel('Force Reaction Y')
    hold on
    plot(armPerturbance(:,2),'g')

    subplot(3,2,5)
    grid on
    ylabel('Force Reaction Z')
    hold on
    plot(armPerturbance(:,3),'b')

    subplot(3,2,2)
    grid on
    ylabel('Torque Reaction X')
    hold on
    plot(armPerturbance(:,4),'r')

    subplot(3,2,4)
    grid on
    ylabel('Torque Reaction Y')
    hold on
    plot(armPerturbance(:,5),'g')

    subplot(3,2,6)
    grid on
    ylabel('Torque Reaction Z')
    hold on
    plot(armPerturbance(:,6),'b')
    

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

    subplot(3,1,2)
    grid on
    ylabel('Torque Reaction Y')
    hold on
    plot(totalTau(:,2),'g')

    subplot(3,1,3)
    grid on
    ylabel('Torque Reaction Z')
    hold on
    plot(totalTau(:,3),'b')

  
    drawnow

end

function graphFval(allData, currSample,pos,numFig)

    fval = zeros(currSample,1);
    for i = 1:currSample
      fval(i) = allData(i).fval;
    end
    
   
    fig = figure(numFig);
    set(fig,'OuterPosition',pos) 
  
    bar(fval(:,1),'b')
    grid on
    drawnow

end

function graphArmStates(allData, currSample,pos,numFig)

    armStates = zeros(currSample,18);
    
    for i = 1:currSample
      armStates(i,:) = allData(i).armStates;
    end
    figure(numFig) 
    
    
        
        for k = 1 : 6

            for j = 1:3

                subplot(6,3,j+(k-1)*3)
                grid on
                hold on
                if j == 1 

                    plot(armStates(:,k),'r') 

                end
                if j == 2

                    plot(armStates(:,k+6),'g')

                end
                if j == 3

                    plot(armStates(:,k+12),'b')

                end

            end

        end 
   
    
    drawnow
end
    
