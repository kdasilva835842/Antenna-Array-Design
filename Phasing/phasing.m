%creating script to generate phasing text file (square array of antennas)
f = 440*10^6;
d = 0.5899;
lam = (3*10^8)/f;
angleDeg = 0:1:180;
phaseDelayDeg =rad2deg(-2*pi*(d/lam).*cos(deg2rad(angleDeg)));
elemPerRow = 10;

phaseDelayVector = zeros(1,elemPerRow^2);
iterations = 0;
for i=1:elemPerRow
    for j=1:elemPerRow
        phaseDelayVector(1,j+iterations*elemPerRow) = phaseDelayDeg(96)*(i-1);
    end
    iterations=iterations+1;
end

magnitude = ones(1,elemPerRow^2);

importVector = [magnitude;phaseDelayVector];
%importVector = transp(importVector);

fileID = fopen('Az180Theta5.txt','w');
fprintf(fileID,'%f %f\r\n',importVector); % '%6.2f %12.8f\r\n'
fclose(fileID);


%%
%RCS vs Altitude relationships

k=1.38064852*10^-23;
T=315; %assuming the system temperature is 45 degrees celcius
B=8026; %system bandwidth in Hz. Check if it must be bandwidth or spectrum analyzer bandwidth
Nf = k*T*B;

syms rAngular positive;

altitude = 1800*10^3;
earthRadius = 6371*10^3;
angle = 0;
earthCenterToAltitude = earthRadius+altitude;
incidentAngle = 90+(90-angle);
eqn = earthCenterToAltitude^2==(rAngular^2)+(earthRadius^2)-2*rAngular*earthRadius*cos(deg2rad(incidentAngle));

distanceToTarget = vpa(solve(eqn,rAngular))/1000

Pt=1000000; %tramsitted power
G=45048.2; %array gain. Transmitter and receiver gain are equal. 38245
%G=33023.6;
r=1800*10^3; %distance from array to target (m)
sphereRadius = 0.05; % radius of spherical debris (m)
RCS=9*pi*((sphereRadius)^2)*(2*pi*sphereRadius/lam)^4; %im assuming rayleigh scattering like the mit paper
Pr=(Pt*(G^2)*RCS*lam^2)/(((4*pi)^3)*r^4);
SNR = 10*log10(Pr/Nf);



%%

noPerRow = [2; 4; 6; 8; 10; 12; 14; 16; 18; 20; 22; 24; 26; 28; 30; 32; 64];
%optimised values for d and h
gainOptVals = [30.246; 137.584; 320.38; 574.068; 904.089; 1312.82; 1794.5; 2348.67; 2980.69; 3689.5; 4470.5; 5325.85; 6259.28; 7267.78; 8348.65; 9505.43; 38245];
% figure(1)
% plot(noPerRow,gain)

optimizedTrend = fit(noPerRow,gainOptVals,'poly2');
figure
plot(optimizedTrend,noPerRow,gainOptVals)

%trends for d=34 and h=17
gainValsSmall = [14.3947; 53.0614; 116.59; 205.177; 318.793; 457.453; 621.153; 809.888; 1023.67; 1262.47; 1526.32; 1815.2; 2129.12; 2468.07; 2832.06; 3221.08];
noPerRowSmall = [2; 4; 6; 8; 10; 12; 14; 16; 18; 20; 22; 24; 26; 28; 30; 32];

smallTrend = fit(noPerRowSmall,gainValsSmall,'poly2');
figure
plot(smallTrend,noPerRowSmall,gainValsSmall)

%trends for d=46.5 and h=15
gainValsMedium = [24.4081; 92.9647; 211.635; 374.549; 584.496; 843.342; 1146.27; 1497.57; 1896.24; 2339.65; 2832.03; 3370.45; 3954.7; 4587.78; 5266.05; 5991.36];
mediumTrend = fit(noPerRowSmall,gainValsMedium,'poly2')
figure
plot(mediumTrend,noPerRowSmall,gainValsMedium)

%plotting the different cases with equations
numRows = 2:1:80;
optCurve = 9.391*(numRows.^2)-3.469*numRows+1.0103;

smallCurve = 3.129*(numRows.^2)+0.4868*numRows+0.9812;

medCurve = 5.853*(numRows.^2)-0.07717*numRows+0.6373;
figure
plot(numRows,0.5*mag2db(optCurve),numRows,0.5*mag2db(medCurve),numRows,0.5*mag2db(smallCurve))
legend('Optimised Case', 'Medium Case', 'Small Case')
grid on;
xlabel('Number of elements per row/column')
ylabel('Gain (dBi)')


