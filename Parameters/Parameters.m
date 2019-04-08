%% Define parameters

%% Fuel Parameters
cFuel1='Gasoline';              % Define fuel types for NASA-database
cFuel2='C2H5OH';

PriceGas = 1.64;                 % Price Gasoline [euro/L]
PriceEth = 1.56;                 % Price Ethanol  [euro/L]

FracGas = 1.0;                  % Mass fraction of gasoline in fuel blend
FracEth = 0.0;                  % Mass fraction of ethanol in fuel blend

%% Model parameters
Throttle = 0.3;
RPM = 3000;               % RPM crankshaft [rotations per minute]
Cw = 2*pi*(RPM/60);             % Angular speed of crankshaft [rad/s]

dCa = 0.005*pi;              % Stepsize crankshaft for calculations [rad]
Numberofcycles = 1;         % Number of times a cycle is completed (one cycle consists of two full rotations)

%% Dimensions engine
B = 68e-3;                % Bore [m]
S = 54e-3;                % Stroke [m]
L = 84e-3;                % Length of connecting rod [m] this might be 84mm
DvC = 8.5;                  % 1 / compression ratio this might be 8.5  
Vd = (S*(1/4)*pi*B^2)/DvC;% Dead volume [m] 

%% Timing parameters
thetaIgn = 2*pi-22.5*pi/180;                % Ignition timing (Start of combustion) [rad]
thetaInOpen = 0*pi;             % Valve timing: Opening intake valves [rad]
thetaInClose = 1.3*pi;            % Valve timing: Closing intake valves [rad]
thetaOutOpen = 2.7*pi;            % Valve timing: Opening exhaust valves [rad]
thetaOutClose = 4.2*pi;           % Valve timing: Closing exhaust valves [rad]

%% Carberettor parameters
FracAreaCarb = 28.28;             % Area air flow / area fuel line [-]
%Throttle = 0.3;                 % Position throttle [0-1] [-]

%% Air parameters
Xair = [0 0 0.21 0 0 0.79];     % Composition air (Mass fractions) [Gasoline, Ethanol, O2, CO2, H2O, N2]

%% Parameters Surroundings
T0 = 293;           % Temperature surroundings [K]
p0 = 1e5;           % Pressure surroundings [Pa]

%% Parameters combustion
thetaComDur = 60*pi/180;       % Duration of combustion [rad]
a = 3;                      % Parameter Wiebe function
m = 4;                      % Parameter Wiebe function

FracGasBurned = 0.6;            % Factor that describes the non-ideal combustion [0-1]

%% Parameters Heat-loss
C1IntExh = 6.18;            % Constant in Woschni equation during intake and exhaust
C1BeforeIgn = 8;         % Constant in Woschni equation before ignition
C1AfterIgn = 3;          % Constant in Woschni equation after ignition
C2IntExh = 0;               % Constant in Woschni equation during intake and exhaust
C2BeforeIgn = 0;            % Constant in Woschni equation before ignition
C2AfterIgn = 23.24e-3;      % Constant in Woschni equation after ignition

Twall = 400;                % Average wall temperature [K]

%% Parameters Mass flow
DpIn = 20e-3;               % Diameter intake valve [m]
DpOut = 20e-3;               % Diameter exhaust valve [m]
Cd = 0.7;                   % Parameter to compensate for not ideal flow typically in range [0.55:0.8]

%% Define constants
global Runiv;
Runiv = 8.314472;       % Gasconstante [J/K mol]

%% Select all species
iSp = myfind({Sp.Name},{cFuel1,cFuel2,'O2','CO2','H2O','N2'});
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];

%% Calculations Cv and Cp for all elements
TR = [1:1:4000];
for i=1:NSp                                                                 % Compute properties for all species for TR temperature range
    Cv(:,i) = CvNasa(TR,SpS(i));
    Cp(:,i) = CpNasa(TR,SpS(i));
    Hi(:,i) = HNasa(TR,SpS(i));
    Ui(:,i) = UNasa(TR,SpS(i));
    Si(:,i) = SNasa(TR,SpS(i));
end
Gamma = Cp./Cv;

%% Calculation compositions
Yair = (Xair.*Mi)/sum(Xair.*Mi);    % Calculation mass fractions air

Yfuel = [FracGas FracEth 0 0 0 0];     % Composition fuel (Mass fractions) [Gasoline, Ethanol, O2, CO2, H2O, N2]
Xfuel = (Yfuel./Mi)/sum(Yfuel./Mi);    % Calculation molar fractions fuel

%% Calculate AF-ratio
for T = TR
    CvFuel(T) = sum(Cv(T,:).*Xfuel);              % Calculation of Cv and Cp of air and fuel for whole temp-range
    CpFuel(T) = sum(Cp(T,:).*Xfuel);
    CvAir(T) = sum(Cv(T,:).*Xair);              
    CpAir(T) = sum(Cp(T,:).*Xair);
end 

AFratio = AFcarb(FracAreaCarb,CpAir(round(T0)),CvAir(round(T0)),CpFuel(round(T0)),CvFuel(round(T0)),p0,T0);     % Determination AF-ratio using carburettor

%% Calculation properties AirFuel mixture:
Ymix= (Yfuel+Yair*AFratio)/(1+AFratio);         % Mass fractions air-fuel mixture
Xmix = (Ymix./Mi)/sum(Ymix./Mi);                  % Molar fractions air-fuel mixture

for T = TR
    CvMix(T) = sum(Cv(T,:).*Xmix);              % Calculation of Cv and Cp of mixture for whole temp-range
    CpMix(T) = sum(Cp(T,:).*Xmix);
    GammaMix(T) = sum(Gamma(T,:).*Xmix);
    HiMix(T) = sum(Hi(T,:).*Xmix);
    UiMix(T) = sum(Ui(T,:).*Xmix);
    SiMix(T) = sum(Si(T,:).*Xmix);
end 

%% Calculation properties Combusted mixture:
% Combustion gasoline: C7.6H13.1 + 10.875 O2 -> 7.6 CO2 + 6.55 H2O
XreacGas = [-1 0 -10.875 7.6 6.55 0];
YreacGas = XreacGas.*Mi/sum(XreacGas.*Mi);
% Combustion ethanol: C2H5OH + 3 O2 -> 2 CO2 + 3 H2O
XreacEth = [0 -1 -3 2 3 0];
YreacEth = XreacEth.*Mi/sum(XreacEth.*Mi);

% Reaction fuel mixture:
Xreacfuel = (Xfuel(1) * XreacGas + Xfuel(2) * XreacEth);
Yreacfuel = Xreacfuel.*Mi/sum(Xreacfuel.*Mi);

% Reaction with enough oxygen:
Xcomb = (Xmix + (Xmix(1)/-Xreacfuel(1))*Xreacfuel) / sum(Xmix + (Xmix(1)/-Xreacfuel(1))*Xreacfuel);
Ycomb = (Xcomb.*Mi)/sum(Xcomb.*Mi);

% If there is not enough oxygen
if Xcomb(3) < 0
    Xcomb = Xmix + (Xmix(3)/-Xreacfuel(3)) * Xreacfuel / sum(Xmix + (Xmix(3)/-Xreacfuel(3)) * Xreacfuel);
    Ycomb = (Xcomb.*Mi)/sum(Xcomb.*Mi);
end

for T = TR
    CvComb(T) = sum(Cv(T,:).*Xcomb);              % Calculation of Cv and Cp of mixture for whole temp-range
    CpComb(T) = sum(Cp(T,:).*Xcomb);
    GammaComb(T) = sum(Gamma(T,:).*Xcomb);
end

%% Initial conditions
Ca(1) = 0;                  % Initial angle of crankshaft [rad] 
P(1) = p0;                  % Initial pressure in cylinder (equal to surroundings) [Pa]
Tcyl(1) = 500;                  % Initial temperature in cylinder (equal to surroundings) [K]
t(1) = 0;                   % Initial time
i = 1;
V(i) = Vcyl(Ca(i),B,S,L,Vd);   % Initial volume in cylinder [m^3]
mresidu = (P(1)/(Tcyl(1)*(CpComb(round(Tcyl(1)))-CvComb(round(Tcyl(1)))))) * Vd;           % Initially the cyclinder is filled with combustion gasses at Pamb and Tamb, so m = rho * Vd

%% Qlhv
%QlhvGasoline = 44.448*10^6;     % Lower Heat Value Gasoline [J/kg]
%QlhvEthanol = 28.865*10^6;      % Lower Heat Value Ethanol [J/kg]
QlhvGasoline = Qlhv(XreacGas(1,:),T0,Mi,Hi);     % Lower Heat Value Gasoline [J/kg]
QlhvEthanol = Qlhv(XreacEth(1,:),T0,Mi,Hi);      % Lower Heat Value Ethanol [J/kg]

Qlhv = QlhvGasoline*FracGas + QlhvEthanol*FracEth;  %% Lower Heat Value Gas mixture [J/kg]


