function [lambdahat] = lamdacombinedfun(neuron,params,Ie_Neurons,Il_Neurons,lambdatype)

% Calculate firing rates frc,fro
[frc,fro,fro2] = fifun2(neuron,params,Ie_Neurons,Il_Neurons,lambdatype);

% lambdatype:
% 1 = MS only
% 2 = opto only
% 3 = MS + Optogenetics (all excitatory - affects all cells)
% indiscriminately)
% 4 = MS + Optogenetics (all silencing - affects all cells indiscriminately
% 5 = MS only, subthreshold
% 6 = opto only, subthreshold
% 7. MS (subthreshold) + Opto (+ Subthreshold) + opto (- supathresh)

% Calculating Lambda hat based off microstimulation
if lambdatype == 1
    lambdahat = frc; % MS only, supathreshold
elseif lambdatype == 2
    lambdahat = fro; % Opto only excitation, supathreshold
elseif lambdatype == 3
    lambdahat = frc + fro; % MS + opto for all, subthreshold
elseif lambdatype == 4
    lambdahat = frc - fro; % MS - opto for all, supathreshold
elseif lambdatype == 5
    lambdahat = frc; % MS
elseif lambdatype == 6
    lambdahat = fro; % MS 
elseif lambdatype == 7
    lambdahat = frc + fro - fro2; % MS + opto for all, subthreshold
end

% Ensuring a firing rate limit is applied
limit = 300;
lambdahat(lambdahat>limit) = limit;
lambdahat(lambdahat<0) = 0;

end