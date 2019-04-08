function Volume = Vcyl(Ca,B,S,L,Vd)                                                  %Function that calculates Volume in cylinder
    CylPosistionMax = 0.5*S + L;                                                    %Position where the cylinder is in its highest position
    CylPosistionCurrent = (0.5*S*cos(Ca) + (L^2 - (0.5*S*sin(Ca)).^2).^0.5);        %Current position of the cylinder
    CylArea = (1/4)*pi*B^2;                                                         %Total area of the cylinder

    Vs = CylArea * (CylPosistionMax - CylPosistionCurrent);                         %The variable volume in the cylinder, so the "not dead" volume

    Volume = Vd + Vs;                                                                    %The total volume in the cylinder
end    