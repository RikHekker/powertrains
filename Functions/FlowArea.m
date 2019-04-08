function A = FlowArea(Ca,thetaOpen,thetaClose,Dp)

D = Dp * (cos((Ca - thetaOpen)/(thetaClose - thetaOpen)*pi))^2;

A = (pi/4) * (Dp^2 - D^2);
