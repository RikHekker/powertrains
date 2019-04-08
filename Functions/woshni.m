function Qdot = woshni (Tcyl,P,Ca,Cv,Cp,C1,C2,Tr,Pr,Vr,Twall,B,S,L,Vd,Cw)
% Calculate heat flow from cylinder into cylinder wall using the woschni equation
gamma = Cp/Cv;  % Calculate gamma

V = Vcyl(Ca,B,S,L,Vd);      % Calculate volume inside cylinder

A = AreaCyl(Ca,B,S,L);      % Calculate contact area cylinder

Pm = Pr * (Vr/V)^gamma;   % Calculate motoring pressure

vp = S * Cw / pi;       % Calculate mean velocity piston

Vs = Vcyl(pi,B,S,L,Vd) - Vd;    % Swept volume is max volume - dead volume

u = C1 * vp + C2 * ((Vs*Tr)/(Pr*Vr)) * (P - Pm);    % Calculate speed in cylinder

h = 3.26 * (P*10^-3)^(0.8) * u^(0.8) * B^(-0.2) * Tcyl^(-0.55);    % Calculate heat transfer rate

Qdot = A * h * (Tcyl - Twall);          % Calculate heat loss to wall