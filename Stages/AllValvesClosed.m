%% During this stage all valves are closed

    for T = TR      % Calculation of matrix of internal energy
    UiCyl(T) = sum(Ui(T,:).*Xcyl);
    end
        
%% Determine reference pressure, temperature and volume for heat flow for
% rest of cylce
Pr = P(i); Tr = Tcyl(i); Vr = V(i);

%% Store composition of gasses, needed for combustion
YcylPreComb = Ycyl;

while Ca(i) <= thetaOutOpen + 4*pi*N                        % Loop for crank angles untill exhaust valve opens
    %% Calculations crank-angle, Volume and time
    i = i + 1;                              % Add one to counter
    waitbar(i / steps)
    Ca(i) = Ca(i-1) + dCa;                  % Add step crankshaft
    V(i) = Vcyl(Ca(i),B,S,L,Vd);            % Calculate new volume
    t(i) = t(i-1) + dCa/Cw;                 % Time is increased by timestep (Stepsize/time per rad)
    MCyl(i) = MCyl(i-1);                    % Mass conservation
    
   %% Calculate mass fraction burned gasses (xb) using the Wiebe function
    if Ca(i)>=thetaIgn + 4*pi*N                              % If the crank angle is greater than the ignition angle combustion takes or has taken place
        xb(i) = FracGasBurned * Wiebe(Ca(i),thetaIgn,thetaComDur,a,m);    % Calculate mass fraction combusted gass using Wiebe function
    else
        xb(i) = 0;                                     % If ignition has not taken place yet no gasses are burned in this cycle
    end
    
    
    %% Calculate composition of gas mixture in cylinder 
    Ycyl = ((1-xb(i)).*(mresidu.*Ycomb + MIntakeThisCycle.*Ymix) + xb(i) .* Ycomb)./sum((1-xb(i)).*(mresidu.*Ycomb + MIntakeThisCycle.*Ymix) + xb(i) .* Ycomb);
    Xcyl = (Ycyl./Mi)./sum(Ycyl./Mi);
    %%

    for T = TR
    CvCyl(T) = sum(Cv(T,:).*Xcyl);              % Calculation of Cv and Cp of mixture in cylinder for whole temp-range
    CpCyl(T) = sum(Cp(T,:).*Xcyl);
    GammaCyl(T) = sum(Gamma(T,:).*Xcyl);
    UiCyl(T) = sum(Ui(T,:).*Xcyl);
    SiCyl(T) = sum(Si(T,:).*Xcyl);
    HiCyl(T) = sum(Hi(T,:).*Xcyl);
    end 
    
    %% Calculate the Temperature and pressure due to combustion
   
     if thetaIgn + 4*pi*N <= Ca(i)<= thetaIgn+thetaComDur + 4*pi*N      
         Tcyl(i)= (Qlhv*(YcylPreComb(1)+YcylPreComb(2))*MCyl(i)*(xb(i)-xb(i-1)))/(CvCyl(round(Tcyl(i-1)))*MCyl(i))+Tcyl(i-1);
         P(i) = (MCyl(i) * (CpCyl(round(Tcyl(i)))-CvCyl(round(Tcyl(i)))) * Tcyl(i))/V(i-1);            % Pressure using gaslaw
     else
       % There are no gasses combusted so the temperature does not change due to combustion
       Tcyl(i) = Tcyl(i-1);
       P(i) = P(i-1);
     end  
 
    %% Calculate effects of compression / expansion
        C1 = P(i)*V(i-1)^(CpCyl(round(Tcyl(i)))/CvCyl(round(Tcyl(i))));         % Calculate Poisson-constant
        P(i) = C1/(V(i)^(CpCyl(round(Tcyl(i)))/CvCyl(round(Tcyl(i)))));       % Calculate pressure using Poisson-relations
        Tcyl(i) = P(i)*V(i)/(MCyl(i) * (CpCyl(round(Tcyl(i)))-CvCyl(round(Tcyl(i)))));           % Calculate temperature using gaslaw
        
        %% Calculate heat flow
        
        % Determine Constants for Woschni equation
        if Ca(i)< thetaIgn + 4*pi*N
            C1Wos = C1BeforeIgn;  C2Wos = C2BeforeIgn;
        else
            C1Wos = C1AfterIgn;  C2Wos = C2AfterIgn;
        end
        % Calculate value heat loss
        Qloss = woshni (Tcyl(i-1),P(i),Ca(i),CvCyl(round(Tcyl(i-1))),CpCyl(round(Tcyl(i-1))),C1Wos,C2Wos,Tr,Pr,Vr,Twall,B,S,L,Vd,Cw) * dCa/Cw;
        
        %% Calculate effects heat loss
        SCyl1 = SiCyl(round(Tcyl(i))).*MCyl(i);      % Determine entropy
        dS = Qloss/Tcyl(i);                 % Determine entropy change due to heat loss
        Tcyl(i) = interp1(SiCyl,TR,(SCyl1+dS)/MCyl(i));                       % Temperature using entropy difference
        P(i) = (MCyl(i) * (CpCyl(round(Tcyl(i)))-CvCyl(round(Tcyl(i)))) * Tcyl(i))/V(i);            % Pressure using gaslaw
      
end