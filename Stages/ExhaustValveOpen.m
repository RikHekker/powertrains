%% During this stage the intake valve is opened

%% Calculation intitial composition gasses in cylinder


for T = TR
    CvCyl(T) = sum(Cv(T,:).*Xcyl);              % Calculation of Cv and Cp of mixture in cylinder for whole temp-range
    CpCyl(T) = sum(Cp(T,:).*Xcyl);
    HiCyl(T) = sum(Hi(T,:).*Xcyl);
    UiCyl(T) = sum(Ui(T,:).*Xcyl);
    SiCyl(T) = sum(Si(T,:).*Xcyl);
end 

while Ca(i) < 4*pi + 4*pi*N                        % Loop for exhaust
    %% Calculations crank-angle, Volume and time
    i = i + 1;                              % Add one to counter
    waitbar(i / steps)
    Ca(i) = Ca(i-1) + dCa;                  % Add step crankshaft
    V(i) = Vcyl(Ca(i),B,S,L,Vd);            % Calculate new volume
    t(i) = t(i-1) + dCa/Cw;                 % Time is increased by timestep (Stepsize/time per rad)
  
    Papprox = (MCyl(i-1)*(CpCyl(round(Tcyl(i-1)))- CvCyl(round(Tcyl(i-1))))*Tcyl(i-1))/(V(i));     % approximation of pressure used to calculate mass flow
    
    %% Calculation mass flow into cylinder
    AValve = FlowArea(Ca(i),thetaOutOpen,thetaOutClose,DpOut);      % Calculate area air can flow through
    
    Mdot(i) = MassflowExhaust(Papprox,p0,Tcyl(i-1),T0,(CpCyl(round(Tcyl(i-1)))/CvCyl(round(Tcyl(i-1)))),AValve,(CpCyl(round(Tcyl(i-1)))- CvCyl(round(Tcyl(i-1)))),Cd);
    MCyl(i) = MCyl(i-1) + Mdot(i) * dCa/Cw;         % Mass in cylinder is increased by mass flow * timestep
     
        
    %% Calculation Heat flow
   Qloss = woshni (Tcyl(i-1),P(i-1),Ca(i-1),CvCyl(round(Tcyl(i-1))),CpCyl(round(Tcyl(i-1))),C1IntExh,C2IntExh,Tr,Pr,Vr,Twall,B,S,L,Vd,Cw)*dCa/Cw;  % Here should the heat flow be implemented, look out Qloss is a heat flow so in J/s!!!!        
    
    %% Calculation Temperature and pressure
    UCyl(i) = UiCyl(round(Tcyl(i-1)))*MCyl(i-1);       % Calculate internal energy at the start of this step
    dU = -Qloss - (P(i-1)-p0)*(V(i)-V(i-1)) + (Mdot(i)*dCa/Cw) * (UiCyl(round(Tcyl(i-1))) + (CpCyl(round(Tcyl(i-1)))-CvCyl(round(Tcyl(i-1)))) * Tcyl(i-1));
    UCyl(i) = UCyl(i) + dU;
    Tcyl(i) = interp1(UiCyl,TR,UCyl(i)/MCyl(i));                            % Calculation temperature using entalphy  
    
    
    % It could happen that the itteration for the mass flow is incorrect due to stepsize, what results in a to low pressure and a temperature that is NaN
    if Tcyl(i) > 0    % If Tcyl is calculated correctly the calculations can proceed          
    P(i) = (MCyl(i) * (CpCyl(round(Tcyl(i)))-CvCyl(round(Tcyl(i)))) * Tcyl(i))/V(i);            % Pressure using gaslaw
    % Calculation of real temperature using Poisson-relations T2 = T1*(P2/P1)^((gamma-1)/gamma)
    %Tcyl(i) = Tcyl(i-1) * (P(i)/P(i-1))^(((CpCyl(round(Tcyl(i-1)))/CvCyl(round(Tcyl(i-1))))-1)/(CpCyl(round(Tcyl(i-1)))/CvCyl(round(Tcyl(i-1)))));
    C1 = Tcyl(i)^(CpCyl(round(Tcyl(i)))/CvCyl(round(Tcyl(i)))) / (P(i-1)^(CpCyl(round(Tcyl(i)))/CvCyl(round(Tcyl(i))) - 1));
    Tcyl(i) = (C1 * P(i)^(CpCyl(round(Tcyl(i)))/CvCyl(round(Tcyl(i))) - 1))^(1/(CpCyl(round(Tcyl(i)))/CvCyl(round(Tcyl(i)))));
    else                % If Tcyl = NaN, the mass flow has to be altered, this is done in the next loop
    P(i) = 0;
    end
    
    if P(i)< p0         % If the pressure drops below surrounding pressure (this is not possible during this stage) the massflow has to be altered.
        
    for factor = 1e-3:0.1:100      % Loop that tries different factors to reduce the mass flow so it does not create a pressure lower than the surrounding pressure
    Papprox = (MCyl(i-1)*(CpCyl(round(Tcyl(i-1)))- CvCyl(round(Tcyl(i-1))))*Tcyl(i-1))/(V(i));     % approximation of pressure used to calculate mass flow
    Papprox = (Papprox - p0)/factor + p0;     % Approximated pressure is changed using the factor defined above
    % Now the same calculations as before this loop are excecuted   
    % Calculation mass flow into cylinder
    AValve = FlowArea(Ca(i),thetaOutOpen,thetaOutClose,DpOut);  % Calculate area air can flow through
    
    Mdot(i) = MassflowExhaust(Papprox,p0,Tcyl(i-1),T0,(CpCyl(round(Tcyl(i-1)))/CvCyl(round(Tcyl(i-1)))),AValve,(CpCyl(round(Tcyl(i-1)))- CvCyl(round(Tcyl(i-1)))),Cd);
    MCyl(i) = MCyl(i-1) + Mdot(i) * dCa/Cw;         % Mass in cylinder is increased by mass flow * timestep
   
    % Calculation Heat flow
    Qloss = 0;
    
    %U(i) = UiCyl(round(Tcyl(i-1)))*MCyl(i-1);
    %dUdt = P(i-1)*((V(i)-V(i-1))/(dCa/Cw)) - Mdot(i)*(-HiCyl(round(Tcyl(i-1))) + Tcyl(i-1)*SiCyl(round(Tcyl(i-1)))) - Qloss;         % Second law of thermodynamics
    %U(i) = U(i) - dUdt * (dCa/Cw);
    %Tcyl(i) = interp1(UiCyl,TR,U(i)/MCyl(i));
    
    UCyl(i) = UiCyl(round(Tcyl(i-1)))*MCyl(i-1);       % Calculate internal energy at the start of this step
    dU = Qloss - (P(i-1)-p0)*(V(i)-V(i-1)) + (Mdot(i)*dCa/Cw) * (UiCyl(round(Tcyl(i-1))) + (CpCyl(round(Tcyl(i-1)))-CvCyl(round(Tcyl(i-1)))) * Tcyl(i-1));
    UCyl(i) = UCyl(i) + dU;
    Tcyl(i) = interp1(UiCyl,TR,UCyl(i)/MCyl(i));                            % Calculation temperature using entalphy
    
    if Tcyl(i) > 0  % Check again if the pressure is good (if NaN this loop has to be runned again)
    P(i) = (MCyl(i) * (CpCyl(round(Tcyl(i)))-CvCyl(round(Tcyl(i)))) * Tcyl(i))/V(i);            % Pressure using gaslaw
    Tcyl(i) = Tcyl(i-1) * (P(i)/P(i-1))^(((CpCyl(round(Tcyl(i-1)))/CvCyl(round(Tcyl(i-1))))-1)/(CpCyl(round(Tcyl(i-1)))/CvCyl(round(Tcyl(i-1)))));
    else
    P(i) = 0;
    end
    if P(i) >= p0    % If the pressure is good this loop can be stopped and the right temperature and pressure are determined.
        break       % Stops only this loop
    end
    end
    end
   
end
%% Determine residue for next cycle
mresidu = MCyl(i);