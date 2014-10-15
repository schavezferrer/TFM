
function sol = tomLab (Ain, bin, lbin, ubin, u0in)

    Name = 'NPSOL SOLUTION';
    
    A       = Ain;          %  Linear constraint
    b_L     = -inf;         %  Lower bound 
    b_U     = bin;          %  Upper bound
    c_L     = [];    %  Two nonlinear inequality constraints
    c_U     = [];           %  Empty means Inf (default)
    x_0     = u0in;         %  Initial value
    x_L     = lbin;         %  Lower bounds on x (signal control u)
    x_U     = ubin;         %  Upper bounds on x (signal control u)
    fLowBnd = 0;            %  A lower bound on the optimal function value
    x_min   = [-2;-2];      %  Used for plotting, lower bounds
    x_max   = [4;4];        %  Used for plotting, upper bounds

    x_opt=[];
    f_opt=[];

    HessPattern = [];  % All elements in Hessian are nonzero. 
    ConsPattern = [];  % All elements in the constraint Jacobian are nonzero. 
    pSepFunc    = [];  % The function f is not defined as separable
    
    
    Prob = conAssign(f, g, H, HessPattern, x_L, x_U, Name, x_0, ...
                                pSepFunc, fLowBnd, ...
                                A, b_L, b_U, c, dc, d2c, ConsPattern, c_L, c_U, ...
                                x_min, x_max, f_opt, x_opt);

%     tomRun('npsol',Prob)

end

% 
%  npsol( A, bl, bu, x, Prob, optPar, ...
%              Warm, H, iState, cLamda, ...
%              SpecsFile, PrintFile, SummFile, ...
%              PriLev, ProbName );