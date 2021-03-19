# C4_dynamic_model
Model description
A dynamic photosynthesis model for C4 crops
This dynamic systems model of C4 photosynthesis was developed based on the previous NADP-ME metabolic model for maize (Wang et al., 2014 ab). The NADP-ME metabolic model is an ordinary differential equation model including all individual steps in C4 photosynthetic carbon metabolism. Here, we extended the model to include posttranslational regulation and temperature response of enzyme activities, dynamic stomata conductance, and leaf energy balance.
This model is written in matlab (R2019a)
Reference 
Wang Y, Chen KX, Long SP. (2021) Toward a Dynamic Photosynthesis Model to Guide the Yield Improvement in C4 Crops. The Plant Journal. 
Wang Y, Bräutigam A, Weber APM, Zhu X-G (2014a) Three distinct biochemical subtypes of C4 photosynthesis? – A modeling analysis. Journal of Experimental Botany 65(13):3567-78.
Wang Y, Long SP, Zhu X-G (2014b) Elements Required for an Efficient NADP-ME Type C4 Photosynthesis-- Exploration using a systems model of C4 photosynthesis. Plant Physiology 164(4):2231-46
Input variables
Photosynthetic parameters
Input file 1: Input_parameter_M.txt
Parameters	Definition	Unit
Slope	Ball-Berry slope	mol m-2 s-1
Intercept	Ball-berry intercept	mol m-2 s-1
ki	Rate constant of stomatal conductance increasing 	min-1
kd 	Rate constant of stomatal conductance decreasing 	min-1
Vpmax	Maximum PEPC activity	µmol m-2 s-1
Vcmax	Maximum rubisco activity	µmol m-2 s-1
taoRubisco	Time constant of Rubisco activation	min
FactorVP	The slope of the linear relationship between measured Vpmax and the maximal PEPC activities in the model	Unitless
FactorVC	The ratio between measured Vcmax and the maximal Rubisco activities in the model	Unitless
RP	PPDK regulatory protein concentration	mmol/L
Rd	Mitochondria respiration	µmol m-2 s-1

Environmental parameters
Input file 2: Input_parameter_Env.txt
Parameters	Definition	Unit
Temperature	Air temperature	°C
CO2_air	CO2 concentration in the air	µmol mol-1
O2_air	O2 concentration in the air	Mmol mol-1
RH	Relative humidity	Unitless
Radiation_NIR	Near-infrared radiation	µmol m-2 s-1
Radiation_LW	Long wave radiation	µmol m-2 s-1
Radiation_PAR	Photosynthetic active radiation	µmol m-2 s-1

Running
Simulating the induction of CO2 uptake rate: Run RAC4leafMetaDrive_Alone.m
Sensitivity analysis: Run Sensitivity.m
Versioning
This is version 1 released on March 2021.
Author
Yu Wang - https://github.com/yuwangcn
