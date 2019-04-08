function A = AreaCyl (Ca,B,S,L)

    CylPosistionMax = 0.5*S + L;                                                    %Position where the cylinder is in its highest position
    CylPosistionCurrent = (0.5*S*cos(Ca) + (L^2 - (0.5*S*sin(Ca)).^2).^0.5);        %Current position of the cylinder
    
    OutlineCyl = B * pi;        % Outline of cylinder is bore * pi
    TopAndBottom = 2 * B^2 * pi/4;  % Surface top and bottom cylinder
    A = OutlineCyl * (CylPosistionMax - CylPosistionCurrent) + TopAndBottom;