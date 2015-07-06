clear

% load computed rates
load data
xe = data(:,1);
phHI = data(:,2);
photogamma = data(:,3);

% set relevant constants
rho = 2.9475329209e-26;
c = 2.99792458e10;
h = 6.6260693e-27;
ev2erg = 1.60217653e-12;
nu0 = 13.6 * ev2erg / h;
F = 1e-2;
intChiE        = 1.0000000000000000e+00;
intChiESigHI   = 1.1391549800556290e-23;
intChiESigHInu = 4.7111661542618619e-41;
TimeUnits = 3.1557e16;
LenUnits = 2.036562e22;
VelUnits = LenUnits / TimeUnits;
DenUnits = 1.0e-26;
aUnits = 1.0;
a = 1.0;
mp = 1.67262171e-24;
dom = DenUnits*a*a*a/mp;
tbase1 = TimeUnits;
xbase1 = LenUnits/a/aUnits;
dbase1 = DenUnits*a*a*a*aUnits*aUnits*aUnits;
coolunit = aUnits*aUnits*aUnits*aUnits*aUnits * xbase1*xbase1 ...
    * mp*mp / tbase1/tbase1/tbase1 / dbase1;
rtunits = ev2erg/TimeUnits/coolunit/dom;


% set derived constants
E = F/c;

% set derived fields
xHI = 1;%ones(size(xe))-xe;
YHI = 0.3908 * (1 - xe.^0.4092).^1.7592;
Ygamma = 0.9971 * (1 - (1 - xe.^0.2663).^1.3163);
phHI_true = c*E*YHI/(h*nu0)*intChiESigHI/intChiE*TimeUnits;
photogamma_true = c*E*xHI.*Ygamma*(intChiESigHI - nu0*intChiESigHInu)/intChiE ...
    * TimeUnits/VelUnits/VelUnits/mp/rtunits;
    
% compute errors
pHI_error = norm((phHI_true - phHI)./phHI_true)
photogamma_error = norm((photogamma_true - photogamma)./photogamma_true)

% plot rates
figure(1)
subplot(2,2,1)
plot(phHI_true,'r-')
xlabel('time step'), ylabel('phHI (analytical)')
subplot(2,2,2)
plot(phHI,'b-')
xlabel('time step'), ylabel('phHI (computed)')

subplot(2,2,3)
plot(photogamma_true,'r-')
xlabel('time step'), ylabel('photogamma (analytical)')
subplot(2,2,4)
plot(photogamma,'b-')
xlabel('time step'), ylabel('photogamma (computed)')

% plot errors
figure(2)
subplot(2,1,1)
plot(phHI_true-phHI)
xlabel('time step'), ylabel('phHI error')

subplot(2,1,2)
plot((photogamma_true-photogamma)./photogamma_true)
xlabel('time step'), ylabel('photogamma error')
