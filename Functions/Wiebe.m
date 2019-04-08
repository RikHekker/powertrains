function xb = Wiebe(ca,thetaIgn,thetaComDur,a,m)
% Wiebe function: combustion model-> calculates mass fraction of burned
% gasses


xb = 1 - exp(-a*((ca-thetaIgn)/thetaComDur)^(m+1));


% Parameters of this equation:
% a: typical 3
% m: typical 4
% theta0: 2*pi-(16/180*pi) [rad]
% deltatheta: 36/180*pi [rad]
end