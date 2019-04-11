%%%%%%%%%%%%%%%%%%%%
% infectivityConstants.m
% Searches for the minimum infectivity s.t. a virus can
% operate to varying degrees of success in a given period of time.
%
% By J. Freedland, 2019.
%
% INPUTS:   params:             Structure with each variable defined.
%                               ex: params.beta == variable beta in the differential equation.
%           percentRange:       Percents at which the virus is deemed complete.
%                               ex: [50 100] considers a virus that is 50, 100% complete
%                               after [timeFrame] hours.
%           timeFrame:          Hours over which percentRange searches.
%
%                               ex: percentComplete = 1:100;
%                                   timeFrame = 12;
%                               The function will search for the infectivity in which 1% -
%                               100% of the virus is spread within 12 hours.
%
% OUTPUTS:  percent:                1xn vector of percents.
%           infectivityConstant:    1xn vector for how many times infective the virus must be.
%%%%%%%%%%%%%%%%%%%

function [percent, infectivityConstant] = infectivityConstants(params,percentRange,timeFrame)

    percent = percentRange;
    infectivityConstant = zeros(1,length(percentRange));
    
    counter = 1;
    for a = percent
        infectivityConstant(1,counter) = findInfectivity(params,a,timeFrame);
        counter = counter+1;
    end

end