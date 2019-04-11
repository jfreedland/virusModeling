%%%%%%%%%%%%%%%%%%%%
% virusTransmission.m
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
%
% OUTPUTS:  t: Vector containing time course. Units depend on rate constants.
%           T: Target cells: Number of healthy cells (as function of t).
%           I: Infected cells: Number of infected cells (as function of t).
%           V: Viruses: Number of viruses (as a function of t).
%%%%%%%%%%%%%%%%%%%

function [t, T, I, V] = virusTransmission(params)

    % Cannot encounter negative concentrations of cells.
    options = odeset('NonNegative', [1 2 3]);

    % Solve differential equation             
    [t1, z1] = ode15s(@infection,[0 params.time_phase],[params.target,params.infected,params.virus,params.lambda,params.dT,params.dI,params.beta,params.p,params.c],...
                                                options);

    % Organize
    T  = (z1(:,1)');
    I = (z1(:,2)');
    V  = (z1(:,3)');
    t = t1';
    
    % Set minimum to 1 cell.
    T(T < 10^0) = 10^0;
    I(I < 10^0) = 10^0;
    V(V < 10^0) = 10^0;
 
    % Models infectivity
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

        % System of differential equations
        dydt(1) = lambda - dT*T - beta*V*T;
        dydt(2) = beta*V*T - dI*I;
        dydt(3) = p*I-c*V;

        % These are the model parameters. Since they do not change with time,
        % their rate of change is 0
        dydt(4) = 0;
        dydt(5) = 0;
        dydt(6) = 0;
        dydt(7) = 0;
        dydt(8) = 0;
        dydt(9) = 0;
        
    end
end
    
    