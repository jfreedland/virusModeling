%%%%%%%%%%%%%%%%%%%%
% modelSystem.m
% Adapted from a differential equation framework courtesy of Soumya Banerjee: www.cs.unm.edu/~soumya.
%
% By J. Freedland, 2019.
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

%           f0:     results of adaptational factors from virusAdaptation.
%
% OUTPUTS:  t: Vector containing time course. Units depend on rate constants.
%           T: Target cells: Number of healthy cells (as function of t).
%           I: Infected cells: Number of infected cells (as function of t).
%           V: Viruses: Number of viruses (as a function of t).
%%%%%%%%%%%%%%%%%%%

function [t, S, I, V, Z] = modelSystem(params)

    % Cannot encounter negative concentrations of cells.
    options = odeset('NonNegative', [1 2 3]);
    
    % Keep initial values stored
    global initial
    initial.v = params.v;
    initial.c = params.c;
    initial.d = params.d;
    initial.n = params.n;

    % Solve differential equation             
    [t1, z1] = ode15s(@infection,[0 params.time_phase],...
        [params.S,params.I,params.V,params.Z,params.mu,params.alpha,params.a,params.beta,params.v, params.k,...
        params.gamma,params.c,params.d,params.delta,params.n],options);

    % Organize
    S   = (z1(:,1)'); % Target cells
    I   = (z1(:,2)'); % Infected cells
    V   = (z1(:,3)'); % Virus
    Z   = (z1(:,4)'); % Immune cells
    t   = t1';
    
    % Set minimum to 1 cell.
    S(S < 10^0) = 10^0;
    I(I < 10^0) = 10^0;
    V(V < 10^0) = 10^0;
    Z(Z < 10^0) = 10^0;
 
    % Models infectivity
    function dydt = infection(t,y)
        
        global type
        global degredationTime
        global healthResponse
        
        % Type 3 or 4 use slowly declining immune systems
        if type >= 3
            deathRate = 24 / degredationTime;   % Hours until fully dead
            factor = -1;                        % Decline entire trajectory
        end
        
        if type == 4
            if (t > (healthResponse / 24))     % Recovery term
                factor = 1;
            end
        end

        dydt = y;

        % Initial conditions
        S  = y(1); % Target cells
        I  = y(2); % Infected cells
        V  = y(3); % Virus
        Z  = y(4); % Immune cells

        % Values of model parameters passed in y
        mu      = y(5); 
        alpha   = y(6); 
        a       = y(7); 
        beta    = y(8);
        v       = y(9); 
        k       = y(10); 
        gamma   = y(11); 
        c       = y(12);
        d       = y(13);
        delta   = y(14);
        n       = y(15);
        
        % Type 2: No immune cell activity
        if type == 2
            v = 0;
            c = 0;
            d = 0;
            n = 0;
        end

        % System of differential equations
        dydt(1) = mu - alpha * S - a * S * V;
        dydt(2) = a * S * V - beta * I - v * I *Z;
        dydt(3) = k * I - gamma * V - a * S * V;
        dydt(4) = n + c * I + d * I * Z - delta * Z;
        
        % These are the model parameters. Since they do not change with time,
        % their rate of change is 0
        dydt(5) = 0;
        dydt(6) = 0;
        dydt(7) = 0;
        dydt(8) = 0;
        dydt(9) = 0;
        dydt(10) = 0;
        dydt(11) = 0;
        dydt(12) = 0;
        dydt(13) = 0;
        dydt(14) = 0;
        dydt(15) = 0;
        
        % Decrease or increase immune cell efficiency, without dropping
        % below zero.
        if type >= 3
            dydt(9)  = factor * initial.v * (deathRate) * double(y(9) > 0);
            dydt(12) = factor * initial.c * (deathRate) * double(y(12) > 0);
            dydt(13) = factor * initial.d * (deathRate) * double(y(13) > 0);
            dydt(15) = factor * initial.n * (deathRate) * double(y(15) > 0);
        end
        
        % Prevent from recovering beyond clincal values.
        if type == 4
            if dydt(9) > initial.v
                dydt(9) = initial.v;
            end
            if dydt(12) > initial.c
                dydt(12) = initial.c;
            end
            if dydt(13) > initial.d
                dydt(13) = initial.d;
            end
            if dydt(15) > initial.n
                dydt(15) = initial.n;
            end
        end
    end
end
    
    
    
    