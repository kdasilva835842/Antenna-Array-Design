%RCS vs Altitude relationships
f = 440*10^6;
d = 0.5899;
lam = (3*10^8)/f;

k=1.38064852*10^-23;
T=315; %assuming the system temperature is 45 degrees celcius
B=8026; %system bandwidth in Hz. Check if it must be bandwidth or spectrum analyzer bandwidth
Nf = k*T*B;

% syms rAngular positive;
% 
% altitude = 1800*10^3;
% %altitude = linspace(200*10^3,1800*10^3,20);
% earthRadius = 6371*10^3;
% angle = 0;
% earthCenterToAltitude = earthRadius+altitude;
% incidentAngle = 90+(90-angle);
% eqn = earthCenterToAltitude^2==(rAngular^2)+(earthRadius^2)-2*rAngular*earthRadius*cos(deg2rad(incidentAngle));
% 
% distanceToTarget = double(vpa(solve(eqn,rAngular)));
% factor = distanceToTarget./altitude;
altitude = linspace(200*10^3,1800*10^3,20);
earthRadius = 6371*10^3;
angle = 0;
earthCenterToAltitude = earthRadius+altitude;
incidentAngle = 90+(90-angle);
a=1;
b=-2*earthRadius*cos(deg2rad(incidentAngle));
c= -((earthCenterToAltitude.^2)-(earthRadius^2));
rAngular= (-b+sqrt((b^2) - 4*a.*c))/(2*a);



Pt=1000000; %tramsitted power
G=45048.2; %array gain. Transmitter and receiver gain are equal. 38245
%G=33023.6;
%r=1800*10^3; %distance from array to target (m)
%r=linspace(200*10^3,1800*10^3,20);
r = rAngular;
%sphereRadius = 0.05; % radius of spherical debris (m)
sphereRadius = linspace(0.025,0.1,4);
bigSphereRadius = 0.5; % this radius should not be less than 0.5
%RCS=9*pi*((sphereRadius)^2)*(2*pi*sphereRadius/lam)^4; %im assuming rayleigh scattering like the mit paper
RCS=9.*pi.*((sphereRadius).^2).*(2.*pi.*sphereRadius./lam).^4;
RCSBig= pi*bigSphereRadius^2;
Pr=(Pt.*(G.^2).*transp(RCS).*lam.^2)./(((4.*pi).^3).*r.^4);
SNR = 10*log10(Pr/Nf);

PrBig = (Pt.*(G.^2).*transp(RCSBig).*lam.^2)./(((4.*pi).^3).*r.^4);
SNRBig = 10*log10(PrBig/Nf);

GsideLobe = 58.8;
anglesideLobe = 2.99;
incidentAngleSideLobe = 90+(90-anglesideLobe);
bSideLobe=-2*earthRadius*cos(deg2rad(incidentAngleSideLobe));
rAngularSideLobe= (-bSideLobe+sqrt((bSideLobe^2) - 4*a.*c))/(2*a);
PrSideLobe=(Pt.*(GsideLobe.^2).*transp(RCS).*lam.^2)./(((4.*pi).^3).*rAngularSideLobe.^4);
SNRSidelobe = 10*log10(PrSideLobe/Nf);
PrBigSideLobe = (Pt.*(GsideLobe.^2).*transp(RCSBig).*lam.^2)./(((4.*pi).^3).*rAngularSideLobe.^4);
SNRBigSideLobe = 10*log10(PrBigSideLobe/Nf);

figure(1)
plot(altitude/1000,SNR(1,:),'b',altitude/1000,SNR(2,:),'r',altitude/1000,SNR(3,:),'k',altitude/1000,SNR(4,:),'m',altitude/1000,SNRBig,'g')
xlabel('Debris Altitude (km)')
ylabel('SNR (dB)')
grid on;
legend('r=2.5 cm','r=5 cm','r=7.5 cm','r=10 cm','r=1 m')

figure(2)
plot(altitude/1000,SNRSidelobe(1,:),'b--',altitude/1000,SNRSidelobe(2,:),'r--',altitude/1000,SNRSidelobe(3,:),'k--',altitude/1000,SNRSidelobe(4,:),'m--',altitude/1000,SNRBigSideLobe,'g--')
xlabel('Debris Altitude (km)')
ylabel('SNR (dB)')
grid on;

figure(3)
plot(altitude/1000,SNR(1,:),'b',altitude/1000,SNR(2,:),'r',altitude/1000,SNR(3,:),'k',altitude/1000,SNR(4,:),'m',altitude/1000,SNRBig,'g')
hold on;
plot(altitude/1000,SNRSidelobe(1,:),'b--',altitude/1000,SNRSidelobe(2,:),'r--',altitude/1000,SNRSidelobe(3,:),'k--',altitude/1000,SNRSidelobe(4,:),'m--',altitude/1000,SNRBigSideLobe,'g--')
hold off;
grid on;
legend('r=2.5 cm','r=5 cm','r=7.5 cm','r=10 cm','r=1 m')
xlabel('Debris Altitude (km)')
ylabel('SNR (dB)')
yticks(-70:10:70)
ylim([-70 70])






