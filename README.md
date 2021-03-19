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

Environmental parameters

Input file 2: Input_parameter_Env.txt


Running
Simulating the induction of CO2 uptake rate: Run RAC4leafMetaDrive_Alone.m
Sensitivity analysis: Run Sensitivity.m


Versioning
This is version 1 released on March 2021.


Author
Yu Wang - https://github.com/yuwangcn
