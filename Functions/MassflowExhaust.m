function Mdot = MassflowExhaust(pcyl,p0,Tcyl,T0,gamma,AValve,Rcyl,Cd)
% This function calculates the mass flow through the valves 

%if p0 <= pcyl          % Intake stage: Gas mixture is sucked into the cylinder
      
    Mdot = -gamma*(p0/pcyl)^(1/gamma)*((2/(gamma-1))*(1-(p0/pcyl)^((gamma-1)/gamma)))^(1/2)*((AValve*pcyl)/(sqrt(gamma*Rcyl*Tcyl)))*Cd;
%end
if p0 == pcyl          % Intake stage: Gas mixture is sucked into the cylinder
    
    %A = (pi/4)*(DpIn^2);       % Calculates the effictive area of the intake valve
    Mdot = 0;
    %Mdot = -gamma*(pcyl/p0)^(1/gamma)*((2/(gamma-1))*(1-(pcyl/p0)^((gamma-1)/gamma)))^(1/2)*((A*p0)/(sqrt(gamma*Ramb*T0)))*Cd;
end

if p0 > pcyl

     Mdot = 0;
end

% Parameters in this formula:
% Dp = diameter valve
% Ds = diameter rod in valve
% Cd = constant that corrects for non ideal flow typical: [0.55:0.8]

end