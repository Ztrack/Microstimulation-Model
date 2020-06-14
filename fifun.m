function [y1,y2] = fifun(Ie_Neurons,Ir_Neurons)

% y1 = Firing rate change due to current
% y2 = Firing rate change due to Optogenetics

%% Current Intensity
% https://link.springer.com/content/pdf/10.1023/A:1008839312043.pdf
% http://www.math.pitt.edu/~bard/pubs/nc100706.pdf

x1 = Ie_Neurons;
A = 60;
Ic = 1;
y1 = A*sqrt(x1);
y1 = real(y1);

%% Luminous Intensity
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6239949/

x2 = Ir_Neurons;
n = 0.82; % Hill coefficient
Imax = 300; % Maximum FR
k = 0.49; % half maximal light sensitivity
y2 = Imax .* ((x2.^n)./((k.^n)+(x2.^n)));

end