function AF_carb = AFcarb(dA, cpair, cvair, cpfuel, cvfuel, p, T)
% Function that calculates the air to fuel ratio in the carburettor;
rho_fuel = p/((cpfuel-cvfuel)*T);
rho_air =  p/((cpair-cvair)*T); 

AF_carb = dA * sqrt(rho_air/rho_fuel);












