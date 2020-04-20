clearvars; clc
load('spreadata.mat')

h = 1; % Angle step
angles = 0:h:360-h; % Range of angles
origin = find(spreadata==max(max(spreadata))); % Origin indices
[y1,x1] = ind2sub([size(spreadata)],origin);  % Origin x,y
dataradii = length(spreadata)-max(x1,y1);
spreadata = spreadata(y1-dataradii:y1+dataradii,x1-dataradii:x1+dataradii); % Center data around origin
k = 1;
for i = 1:4

    if i == 1
        for ii = 1:length(spreadata)
            indata(k) = sub2ind(size(spreadata),1,ii);
            k = k+1;
        end
    elseif i == 2
        for ii = 1:length(spreadata)
            indata(k) = sub2ind(size(spreadata),length(spreadata),ii);
            k = k+1;
        end
    elseif i == 3
        for ii = 1:length(spreadata)
            indata(k) = sub2ind(size(spreadata),ii,1);
            k = k+1;
        end
    elseif i == 4
        for ii = 1:length(spreadata)
            indata(k) = sub2ind(size(spreadata),ii,length(spreadata));
            k = k+1;
        end
    end
    
end

for i = 1:length(indata)
    Y = []; X = [];
    [y2,x2] = ind2sub([size(spreadata)],indata(i));
    theta = atan2d(y2-y1,x2-x1); % Angle from -180 to 180
    normDeg(i) = mod(theta,360);
end

for i = 1:length(angles)
    closestangle = abs(angles(i)-normDeg);
    bestind(i) = find(closestangle == min(closestangle));
end

for i = 1:length(bestind)
    Y = []; X = [];
    [y2,x2] = ind2sub([size(spreadata)],indata(bestind(i)));
    xc = [x1 x2]; % X line
    yc = [y1 y2]; % Y line
    Points = max(abs(diff(xc)), abs(diff(yc)))+1; % Number of points in line
    rIndex = round(linspace(yc(1), yc(2), Points)); % Row indices
    cIndex = round(linspace(xc(1), xc(2), Points)); % Column indices
    lindex = sub2ind([size(spreadata)], rIndex,cIndex); % Linear indices
    
    for ii = 1:length(lindex)
        [y3,x3] = ind2sub([size(spreadata)],lindex(ii));
        X(ii) = sqrt((x3-x1)^2+(y3-y1)^2); % Distance away from origin for each point
    end
    X(X == 0) = 1;
    Y = spreadata(lindex); % Data value at this point/index
    
    f = fit(X',Y','Power2');
    a(i) = f.a; % 
    b(i) = f.b; % 
    c(i) = f.c; % 
    %d(i) = f.d; %

end

% Plotting
figure;
x = 1:500; % Distance
Y1 = []; Y2 = [];
for i = 1:length(angles)
    Y = a(i).*x.^b(i)+c(i);
    Y(Y<0) = 0;
    plot(Y);
    hold on;
    Y1(i,:) = Y;
end
title('Light Spread Function Fitting');
xlabel('Distance from origin um');
ylabel('Magnitude of Power');

figure;
cangles = 45;
for i = 1:360/45
    % a*x^b+c format
    Y2 = mean(Y1((i-1)*cangles+1:i*cangles,:));
    plot(Y2);
    hold on;
end
title('Light Spread Function Fitting');
xlabel('Distance from origin um');
ylabel('Magnitude of Power');
legend('1-45','46-90','91-135','136-180','181-225','226-270','271-315','316-360');

Y3 = mean(Y1);
lightspread.averaged = fit(x',Y3','Power2');
lightspread.equation = ['a*x^b+c'];
lightspread.all.a = a;
lightspread.all.b = b;
lightspread.all.c = c;

% General Averaged Case
figure; 
subplot(2,1,1); plot(x,Y3);
title('1mW Light Intensity Distribution From Center');
xlabel('Distance From Origin (um)');
ylabel('Relative Light Intensity');
Stim_Loc = zeros(500,500);
Stim_Loc(250,250) = 1;
Ed = bwdist(Stim_Loc); Ed(250,250) = 1;
Stim_int = lightspread.averaged.a.*Ed.^lightspread.averaged.b+lightspread.averaged.c;
subplot(2,1,2); imagesc(Stim_int); colorbar; xlabel('um'); ylabel('um'); title('1mW Light Intensity Distribution');

% Angle Specific
figure;
subplot(2,1,1);
cangles = 45;
for i = 1:360/45
    % a*x^b+c format
    Y2 = mean(Y1((i-1)*cangles+1:i*cangles,:));
    plot(Y2);
    hold on;
end
title('Averaged 1mW Light Intensity Distribution From Center');
xlabel('Distance from origin um');
ylabel('Magnitude of Power');
legend('1-45','46-90','91-135','136-180','181-225','226-270','271-315','316-360');
subplot(2,1,2);
Stim_Loc = zeros(500,500);
Stim_Loc(250,250) = 1;
Ed = bwdist(Stim_Loc); Ed(250,250) = 1;
indata = size(Ed,1)*size(Ed,2); % Indices
Stim_int = zeros(500,500);
for i = 1:length(indata)
    [y2,x2] = ind2sub([size(Ed)],indata(i));
    theta = atan2d(y2-y1,x2-x1); % Angle from -180 to 180
    theta = round(mod(theta,360)); % Normalize angle & round to nearest
    Stim_int(indata(i)) = lightspread.all.a(i).*Ed(indata(i)).^lightspread.all.b(i)+lightspread.all.c(i);
end
Stim_int = lightspread.averaged.a.*Ed.^lightspread.averaged.b+lightspread.averaged.c;
subplot(2,1,2); imagesc(Stim_int); colorbar; xlabel('um'); ylabel('um'); title('1mW Light Intensity Distribution');

% Electrical Stim Spread to compare
figure;
subplot(2,1,1);
X = 1:50; % um
Y = 100./X.^2;
plot(X,Y);
title('Current Intensity Distribution From Center');
xlabel('Distance from origin um');
ylabel('Magnitude of Current (AU)');
subplot(2,1,2);
Stim_Loc = zeros(50,50);
Stim_Loc(25,25) = 1;
Ed = bwdist(Stim_Loc); Ed(25,25) = 1;
Stim_int = zeros(50,50);
Stim_int = 1000./Ed;
subplot(2,1,2); imagesc(Stim_int); colorbar; xlabel('um'); ylabel('um'); title('Current Intensity Distribution');



%save('lightspread.mat','lightspread');