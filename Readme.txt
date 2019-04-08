Readme complex model

The complex model consists of several different files devided in different folders.
Below is an overview of the different files and folders given:

"Functions": This folder contains all the defined used functions that are in a separe m-file
	- AFcarb.m:		Calculates the air-fuel ratio for the given fuel composition
	- AreaCyl.m:		Calculates the contact area of the cylinder with the gasses inside
	- FlowArea.m:		Calculates the flow area for gasses through the intake/exhaust valve
	- MassflowExhaust.m:	Calculates the mass flow per step through the exhaust valve
	- MassflowIntak.m:	Calculates the mass flow per step through the intake valve
	- Qlhv.m:		Calculates the lower heating value for the given fuel composition
	- SpecificHeat.m:	Calculates the specific heat of a mixture
	- Vcyl.m:		Calculates the volume inside the cylinder for a given crank angle
	- Wiebe.m:		The Wiebe function that models the combustion rate
	- woschni.m:		Woschni function that models heat transfer between gas and cylinder wall

"NASA_database": This folder contains the given NASA-polynomals to determine thermodynamic properties of the substances

"Parameters": This folder constains all files that define the parameters
	- Parameters.m:		In this script all the parameters are defined, these can be altered by the user

"Stages": This folder contains the scripts that model the different stages of the thermodynamic cycle
	- AllValvesClosed.m:	This models the phases: compression, expansion and combustion (The stages where all valves are closed)
	- ExhaustValveOpen.m:	This models the exhaust phase (The exhaust valve is opened)
	- IntakeValveOpen.m:	This model the intake phase (The intake valve is opened)

Run_complex_model.m: This is the main script of the model, this runs all other scripts in the right order.

--------------------------------------------------------------------------------------------------------
HOW TO USE THE MODEL:

The model has been build in such way that the user only needs to use two scripts:
- Run_complex_model.m
- Parameters/Parameters.m

All the parameters of the model can be entered in the script Parameters.m.
Then the model can be runned by running the script Run_complex_model.m