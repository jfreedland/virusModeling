%% HIV model
% source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5072357/

% Initial concentrations of cells.
params.target      = 5e3;     % Healthy cells per uL
params.infected    = 0;      % Infected cells per uL
params.virus       = 0.4e-3;  % Viruses per uL

% Virus-specific spread parameters.
params.lambda  = 100;                   % Production of new healthy cells (per uL)
params.dT      = 0.1;                   % Healthy cell death rate (without virus)
params.dI      = 0.5;                   % Infected cell death rate
params.beta    = 10^-5;                 % Virus infectivity (per uL^-1)
params.p       = 1.5*10^3;              % Virus production rate
params.c       = 10;                    % Virus clearence rate (by immune system)
params.omega   = 10^-3;

params.time_phase   = 30;      % Number of hours/days to consider (depends on units)

fSolution = zeros(100,3);
for a = 1:100
    [f0, ~] = virusAdaptation(params,1,1);
    fSolution(a,:) = f0;
    disp(f0)
end
% f0 = [1 1 1];
[t, ~, I, ~] = modelSystem(params,f0);
plot(t,I)
A = [t',I'];
%% Dengue Virus
% source: https://www.sciencedirect.com/science/article/pii/S0895717708002732

% Amount of solution considered (in ml)
mL = 27;

% Initial concentrations of cells.
params.target      = 5e3;     % Healthy cells per uL
params.infected    = 0;      % Infected cells per uL
params.virus       = 0.4e-3;  % Viruses per uL

% Values of model parameters passed in y
params.lambda  = 80;         % per uL
params.dT      = 0;
params.dI      = 0.5; 
params.beta    = 10^-3;  % per uL^-1
params.p       = 20;
params.c       = 0.8;
params.omega   = 10^-3;

params.time_phase   = 30;

% [f0, ~] = virusAdaptation(params,1,1);
f0 = [1 1 1];
[t, ~, I, ~] = modelSystem(params,f0);

% fSolution = zeros(100,3);
% for a = 1:100
%     [f0, ~] = virusAdaptation(params,1,1);
%     fSolution(a,:) = f0;
%     disp(f0)
% end
% success = zeros(100,1);
% A = sort(fSolution);
% 
% for a = 1:100
%     f0 = A(a,:);
%     [t, ~, I, ~] = modelSystem(params,f0);
%     
%     [~,i] = min(abs(t - 1));
%     success(a,1) = I(i) / max(I(:));
% end