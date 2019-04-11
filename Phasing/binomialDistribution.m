elementsPerRow = 76;%***************************
f = 440*10^6;
d = 0.5899;
lam = (3*10^8)/f;
%angleDeg = 0:1:180;
desiredTheta = 5;%*****************************
desiredPhi = 45;%*******************************
%phaseDelayDeg =rad2deg(-2*pi*(d/lam).*cos(deg2rad(angleDeg)));
phaseDelayY =rad2deg(-2*pi*(d/lam)*sin(deg2rad(desiredTheta))*sin(deg2rad(desiredPhi)));
phaseDelayX =rad2deg(-2*pi*(d/lam)*sin(deg2rad(desiredTheta))*cos(deg2rad(desiredPhi)));
%desiredAngleAz0to180 = 20;%***************************

phaseDelayVectorY = zeros(1,elementsPerRow^2);
phaseDelayVectorX = zeros(1,elementsPerRow^2);
iterations = 0;
for i=1:elementsPerRow
    for j=1:elementsPerRow
        %phaseDelayVector(1,j+iterations*elementsPerRow) = phaseDelayDeg(-desiredAngleAz0to180+91)*(i-1);
        phaseDelayVectorY(1,j+iterations*elementsPerRow) = phaseDelayY*(i-1);
    end
    iterations=iterations+1;
end

iterations = 0;
for i=1:elementsPerRow
    for j=1:elementsPerRow
        phaseDelayVectorX(1,j+iterations*elementsPerRow) = phaseDelayX*(j-1);
    end
    iterations=iterations+1;
end

topLeftQuarter = pascal(elementsPerRow/2);
bottomLeftQuarter = flip(topLeftQuarter);
leftHalf = vertcat(topLeftQuarter,bottomLeftQuarter);
rightHalf = flip(leftHalf,2);
biDistribution = horzcat(leftHalf,rightHalf);

%h = fspecial('gaussian',1);
Y = imgaussfilt(biDistribution,elementsPerRow*0.3241,'FilterDomain','spatial'); %2*ceil(2*SIGMA)+1
%Y = filter2(h,biDistribution);
Y = Y./max(max(Y));
Yrounded = round(Y,2);

figure(1)
surf(biDistribution)

figure(2)
surf(Y)


%singleColumn = reshape(biDistribution,[1,elementsPerRow^2]);
singleColumn = reshape(Y,[1,elementsPerRow^2]);

importVector = [singleColumn;phaseDelayVectorX + phaseDelayVectorY];

fileID = fopen('binomial.txt','w');
fprintf(fileID,'%f %f\r\n',importVector); % '%6.2f %12.8f\r\n'
fclose(fileID);
%%
maxPowerForElement = 350;
elementPowerMatrix = round(Y.*maxPowerForElement);
totalPowerTransmitted = sum(sum(elementPowerMatrix));

figure
surf(elementPowerMatrix)



