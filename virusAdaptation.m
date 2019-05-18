%%%%%%%%%%%%%%%%%%%%
% virusAdaptation.m
%
% This function solves a series of ODEs for how a virus spreads throughout the body.
% It then evolves three parameters to find parameters for achieving
% a certain percentage of spread within a specific period of time.
%
% By J. Freedland, 2019.
% Uses an ODE framework courtesy of Soumya Banerjee: www.cs.unm.edu/~soumya.
%
% INPUTS:   params: Structure with each variable defined.
%                   ex: params.beta == variable beta (infectivity) in the differential equation.
%       
%       REQUIRES:   params.target       Healthy cells (initial condition)
%                   params.infected     Infected cells (initial condition)
%                   params.virus        Viruses (initial condition)
%                   params.lambda       Production of new healthy cells
%                   params.dT           Healthy cell death rate without virus
%                   params.dI           Infected cell death rate
%                   params.beta         Virus infectivity
%                   params.p            Virus production rate
%                   params.c            Virus clearence rate by immune system
%                   params.time_phase	Number of hours/days to consider (depends on units)
%                   params.omega        Cell-to-cell infectivity
%
%           t1:             amount of time for our virus to successfully
%                           infect percentSpread of the body
%           percentSpread:  percent of infected cells considered successful (0 - 1)
%
% OUTPUTS:  f0:     Adaptation parameters in the form [param.beta param.p params.omega].
%                       ex: if f0 = [1 1 2], doubling the cell-to-cell
%                       infectivity of a virus would allow it to spread.
%           realismFactor: structure with data from fminsearch.
%%%%%%%%%%%%%%%%%%%

function f0 = virusAdaptation(params,t1,percentSpread)

    global trajectory

    f0 = [1 1 1]; % Begin search assuming no virus adaptation.
    
    % Calculate new parameters.
    adaptVirusInfectivity   = @(f) f(1) * params.beta;
    adaptVirusReplication   = @(f) f(2) * params.p;
    adaptCellInfectivity    = @(f) f(3) * params.omega;
    
    % Define differential equations
    options = odeset('NonNegative', [1 2 3]); % Cannot encounter negative concentrations of cells.
    trajectory = @(f) ode15s(@infection,[0 params.time_phase],[params.target,params.infected,params.virus,params.lambda,...
        params.dT,params.dI,adaptVirusInfectivity(f),adaptVirusReplication(f),params.c, adaptCellInfectivity(f)],...
        options);
    
    % We pass our variable to a function to calculate a numerical ratio.
    calculateRatio = @(f) calculateRatioFunc(f,percentSpread,t1);
    minimizeRatio = @(f) abs(1 - calculateRatio(f));     
    
    % Use gradient descent to find minimum
    for iter = 1:20
        f0 = abs(fminsearch(@(f) minimizeRatio(f), f0));
    end
      
    % Differential equations.
    function dydt = infection(t,y)

        dydt = y;

        % Initial conditions
        T  = y(1); % Target cells
        I  = y(2); % Infected cells
        V  = y(3); % Virus

        % Values of model parameters passed in y
        lambda  = y(4); 
        dT      = y(5); 
        dI      = y(6); 
        beta    = y(7);
        p       = y(8); 
        c       = y(9); 
        omega   = y(10); 

        % System of differential equations
        dydt(1) = lambda - omega * I * T - dT*T - beta*V*T;
        dydt(2) = beta*V*T + omega * I * T - dI*I;
        dydt(3) = p*I-c*V;

        % These are the model parameters. Since they do not change with time,
        % their rate of change is 0
        dydt(4) = 0;
        dydt(5) = 0;
        dydt(6) = 0;
        dydt(7) = 0;
        dydt(8) = 0;
        dydt(9) = 0;
        dydt(10) = 0;
    end

    function ratio = calculateRatioFunc(f,percentSpread,t1)  
        
        A       = trajectory(f);          % solve differential equation for variable
        [~,i]   = min(abs(A.x-t1));       % find initial time index

        infectedCells   = A.y(2,i);       % number of infected cells in system
        bestCase        = max(A.y(2,:));  % maximum number of cells infected

        % Calculate percent of cells that have been infected after __ time.
        ratio   = infectedCells / (bestCase .* percentSpread);
    end
end
    
    