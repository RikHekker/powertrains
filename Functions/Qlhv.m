function Qlhv = Qlhv(Xreac, T,Mi,Hi)
%Xreac: the molare ratio between the elements of the reaction
%T: the corresponding temperature

mreac(1,:) = Xreac.*Mi;                          %Calculate the mass
mreacperfuel = mreac(1,:)./ (mreac(1,1)+mreac(1,2));  %Calculate the mass per fuel
Qlhv = sum(mreacperfuel.*Hi(T,:));                    %Calculate the Qlhv




%%Qlhv for pure gasoline:
  %Xreac= XreacGas
  %T=T0

%%Qlhv for pure ethanol:
  %Xreac: XreacEth
  %T=T0