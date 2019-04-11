%%%%%%%%%%%%%%%%%%%%
% findInfectivity.m
% Searches for the minimum infectivity using a virus model s.t. a virus can
% operate in a given period of time.
%
% By J. Freedland, 2019.
%
% INPUTS:   params:             Structure with each variable defined.
%               ex: params.beta == variable beta in the differential equation.
%           percentComplete:    Percent at which we consider the reaction fully completed. 
%           timeFrame:          Hours at which percentComplete searches.
%
%               ex: percentCompelete = 85;
%                   timeFrame = 8;
%               The function will search for the infectivity in which 85%
%               of the virus spreads by 8 hours.
%
% OUTPUTS:  infectivityConstant:    how many times more infective our new
%                                   virus is.
%%%%%%%%%%%%%%%%%%%

function infectivityConstant = findInfectivity(params,percentComplete,timeFrame)

    time = timeFrame / 24; % convert to days
    percent = percentComplete / 100; % convert to decimal
    stableBeta = params.beta;

    search = 0:100:10000;
    
    for a = search % all possible infectivities
        params.beta    = stableBeta*a; % try infectivity
        [t1, T, I, V] = virusTransmission(params);
        [~,i] = min(abs(t1-time)); % find time mark
        
        I(I > max(T)) = max(T); % point at which rxn is ~100% complete

        if V(i) > 10 % must produce virus (successful infection)
            if abs(log10(I(i)) - log10(max(I))) <= log10(max(I)) * (1-percent) % reaches percent

                for b = a - 100:a % 1 resolution
                    params.beta    = stableBeta*b;
                    [t1, T, I, V] = virusTransmission(params);
                    [~,i] = min(abs(t1-time));
                    
                    I(I > max(T)) = max(T);
                    
                    if V(i) > 10
                        if abs(log10(I(i)) - log10(max(I))) <= log10(max(I)) * (1-percent)
                            for c = b - 1:0.01:b % 0.01 resolution
                                params.beta   = stableBeta*c;
                                [t1, T, I, V] =  virusTransmission(params);
                                [~,i] = min(abs(t1-time));
                                
                                I(I > max(T)) = max(T);
                    
                                if V(i) > 10
                                    if abs(log10(I(i)) - log10(max(I))) <= log10(max(I)) * (1-percent)
                                        infectivityConstant = c;
                                        disp(c)
                                        return
                                    end
                                end
                            end
                        end
                    end
                end           
            end
        end
    end
end