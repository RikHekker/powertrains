function [cv, cp, gamma] = SpecificHeat(T,FuelType)
%% Load Nasadatabase
TdataBase=fullfile('NASA_database','NasaThermalDatabase');
load(TdataBase);

%% Nasa is ready
iSp = myfind({Sp.Name},{FuelType});                            % Proper indices to database
cv = CvNasa(T,Sp(iSp));
cp = CpNasa(T,Sp(iSp));
gamma = cp./cv;
end

