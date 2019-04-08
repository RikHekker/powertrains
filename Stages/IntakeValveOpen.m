%% During this stage the intake valve is opened
MIntakeThisCycle = 0;
MCyl(i) = mresidu;

% Calculation intitial composition gasses in cylinder
Xcyl = Xcomb;
Ycyl = (Xcyl.*Mi)./sum(Xcyl.*Mi);

for T = TR
    CvCyl(T) = sum(Cv(T,:).*Xcyl);              % Calculation of Cv and Cp of mixture in cylinder for whole temp-range
    CpCyl(T) = sum(Cp(T,:).*Xcyl);
    GammaCyl(T) = sum(Gamma(T,:).*Xcyl);
    HiCyl(T) = sum(Hi(T,:).*Xcyl);
    UiCyl(T) = sum(Ui(T,:).*Xcyl);
end 
HCyl(i) = MCyl(i) * HiCyl(round(Tcyl(i))) + (P(i)-p0)*(V(i)-V(i));
UCyl(i) = MCyl(i) * UiCyl(round(Tcyl(i)));

% Determine reference pressure, temperature and volume for heat flow
Pr = P(i); Tr = Tcyl(i); Vr = V(i);

while Ca(i) < thetaInClose + 4*pi*N                        % Loop for intake
    % Calculations crank-angle, Volume and time
    i = i + 1;                              % Add one to counter
    waitbar(i / steps)
    Ca(i) = Ca(i-1) + dCa;                  % Add step crankshaft
    V(i) = Vcyl(Ca(i),B,S,L,Vd);            % Calculate new volume
    t(i) = t(i-1) + dCa/Cw;                 % Time is increased by timestep (Stepsize/time per rad)
  
    P(i) = (MCyl(i-1)*(CpCyl(round(Tcyl(i-1)))- CvCyl(round(Tcyl(i-1))))*Tcyl(i-1))/(V(i));
    
    % Conditions in carburettor
    Tcarb = T0;
    
    % Calculation mass flow into cylinder
    AValve = FlowArea(Ca(i),thetaInOpen,thetaInClose,DpIn);    % Calculate area where air can flow through
        
    Mdot(i) = MassflowIntake(P(i),p0,Tcarb,(CpMix(round(Tcarb))/CvMix(round(Tcarb))),AValve,(CpMix(round(Tcarb))- CvMix(round(Tcarb))),Cd)*Throttle;
    MCyl(i) = MCyl(i-1) + Mdot(i) * dCa/Cw;         % Mass in cylinder is increased by mass flow * timestep
    MIntakeThisCycle = MIntakeThisCycle + Mdot(i) * dCa/Cw;    % Total mass intake this cycle
    
    % Calculation composition gasses in cylinder
    Ycyl = (mresidu*Ycomb + MIntakeThisCycle*Ymix)/sum(mresidu*Ycomb + MIntakeThisCycle*Ymix);
    %Xcyl = mfracfuel(MIntakeThisCycle,mresidu,AFratio) * Xfuel + mfracair(MIntakeThisCycle,mresidu,AFratio) * Xair + mfracresidu(MIntakeThisCycle,mresidu) * Xcomb;
    Xcyl = (Ycyl./Mi)./sum(Ycyl./Mi);
    
    for T = TR
    CvCyl(T) = sum(Cv(T,:).*Xcyl);              % Calculation of Cv and Cp of mixture for whole temp-range
    CpCyl(T) = sum(Cp(T,:).*Xcyl);
    GammaCyl(T) = sum(Gamma(T,:).*Xcyl);
    HiCyl(T) = sum(Hi(T,:).*Xcyl);
    UiCyl(T) = sum(Ui(T,:).*Xcyl);
    end 
    
    % Calculation Heat flow
    Qloss = woshni (Tcyl(i-1),P(i-1),Ca(i-1),CvCyl(round(Tcyl(i-1))),CpCyl(round(Tcyl(i-1))),C1IntExh,C2IntExh,Tr,Pr,Vr,Twall,B,S,L,Vd,Cw) * dCa/Cw;
    
    
    %% Calculation Temperature and pressure
    dU = -Qloss + (P(i-1)-p0)*(V(i)-V(i-1)) + (Mdot(i)*dCa/Cw) * (UiMix(round(Tcarb)) + (CpMix(round(Tcarb))-CvMix(round(Tcarb))) * Tcarb);     % First law for open systems
    UCyl(i) = UCyl(i-1) + dU;       % New internal energy
    Tcyl(i) = interp1(UiCyl,TR,UCyl(i)/MCyl(i));                            % Calculation temperature using entalphy
        
    P(i) = MCyl(i) * (CpMix(round(T0))-CvMix(round(T0))) * Tcyl(i) / V(i);  % Calculate pressure using gas law
end