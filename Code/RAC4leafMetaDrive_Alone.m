SpeciesCol=3; %1Maize	2sorghum	3sugarcane
%Enzyme activation settings
global EAPPDK;
global PRac 
global RedoxEnyAct;
global GsResponse;
EAPPDK=1;%%if EAPPDK=0 PPDK is fully actived; %%if EAPPDK=1 include PPDK activation by PDRP
PRac=1;%%if PRac=0 Rubisco is fully actived; %%if PRac=1 include Rubisco activation by Rca
RedoxEnyAct=1; %%if RedoxEnyAct=0 activities of photosynthetic enzymes are not regulated by light intensity; %%if RedoxEnyAct=1 include light induced enzyme activation 
GsResponse=1; %%if GsResponse=0 Ball berry model no time dependent Gs response ; %%if GsResponse=1 include Gs response, using ki and kd

%Inpute parameter settings
global Vexsen;
global parameterdata;
parametertable=importdata('../input/Input_parameter_M.txt');
parameterdata=parametertable.data;
global parameter_Env;
parametertable_Env=importdata('../input/Input_parameter_Env.txt');
parameter_Env=parametertable_Env.data;

global Vpmax;
global Vcmax;
global ki;
global kd;
global taoRub;
global FactorVP;
global FactorVC;
global PPDKRP;
global MRd;
global MeasuredTemp;
ki=parameterdata(3,SpeciesCol);
kd=parameterdata(4,SpeciesCol);
slop=parameterdata(1,SpeciesCol);
inter=parameterdata(2,SpeciesCol);
taoRub=parameterdata(7,SpeciesCol);
FactorVP=parameterdata(8,SpeciesCol);
FactorVC=parameterdata(9,SpeciesCol);
Vpmax=parameterdata(5,SpeciesCol);
Vcmax=parameterdata(6,SpeciesCol);
PPDKRP=parameterdata(10,SpeciesCol);
MRd=abs(parameterdata(11,SpeciesCol));
MeasuredTemp=parameterdata(12,SpeciesCol);
global Tao_MDH;
global Tao_PEPC;
global kdcon;
Tao_MDH=1;
Tao_PEPC=2;
kdcon=1;%=1: constant kd; =0: kd change with light

%Sensitivity parameter settings
global VPPDKsen;
global VMDHsen;
global VMEsen;
global VGAPDHsen;
global VSBPsen;
global VFBPsen;
global VPRKsen;
global Jmaxsen;
global gsxsen;
global Jsen;
VPPDKsen=1;
VMDHsen=1;
VMEsen=1;
VGAPDHsen=1;
VSBPsen=1;
VFBPsen=1;
VPRKsen=1;
Jmaxsen=1;
Jsen=1;
Vexsen=1;
gsxsen=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Relative humidity setting
%Saturated water vapor: ESatweather = 0.611 * exp(17.502 * WeatherTemperature / (WeatherTemperature + 240.97));
%RH=0.65;% Original setting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Para_mata;
global VmaxC4;
global BallBerryInterceptC4;
global BallBerrySlopeC4;
Para_mata=1;%%if Para_mata=1, C4 Metabolic model and Gs model integrated  if Para_mata=0 Steady-state mdoel and gs model
VmaxC4=160;%For steady-state photosynthesis mdoel
BallBerryInterceptC4=inter*1.6;
BallBerrySlopeC4=slop;%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Environmental parameters
global WeatherTemperature;
global Air_CO2;
global Air_O2;
global WeatherRH;
global WeatherWind;
global Radiation_NIR;
global Radiation_LW;
global Radiation_PARo;
global PhiLeaf

%WeatherTemperature=25;
% Air_CO2=410;
% Air_O2=210.0;
% WeatherRH=RH;
% Radiation_NIR=0;
% Radiation_LW=0;
WeatherWind=3.5;%m/s %Set as 3.5 to match the bundary layer conductance from licor data
PhiLeaf=0;%Mpa
WeatherTemperature=parameter_Env(1);
Air_CO2=parameter_Env(2);
Air_O2=parameter_Env(3);
WeatherRH=parameter_Env(4);
Radiation_NIR=parameter_Env(5);
Radiation_LW=parameter_Env(6);
Radiation_PARo=parameter_Env(7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global CI;
global O2;
CI=Air_CO2*0.4/(3 * 10^4);%The initial intercelluar CO2 Concnetration mmol/L
O2=Air_O2*1.26/1000;%O2 concentration 
global phi;
global Lpd;
phi=0.03;
Lpd=400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%C4 decarboxylation pathway settings
global C13ratio;% C insotope simulation
C13ratio=0;%0.01115;
global U;
global V;
U=0;% light partition coefficient
V=0;% light partition coefficient
global Ratio;
Ratio=4; % Enezyme activity variation factor
global pathway_option;
global RatioPPDK;%PPDK in BSC
RatioPPDK=0;
global Pvpr8M;
Pvpr8M=0;%if 0 Glycerate kinase in BSchl; if 1 Glycerate kinase in Mchl;
global Bchl_CP;%Pi content in BS chloroplast
Bchl_CP= 25.0;
global MC_CP;%Pi content in MC
MC_CP=15.0;
global Mchl_CP;%Pi content in MC chloroplast
Mchl_CP=15.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathway_option=0;
%%% 0=Normol NADP-ME type 
%%% 1=Asp+Mal transport and MDH type 
%%% 2=Asp+Mal and PCK type 
%%% 3 Asp+Mal and PCK+MDH type 
%%% 4 Asp and PCK only type
%%% 6 DiT2 mutant
%%% 7 NAD-ME type
%%% 8 NAD-ME+PCK type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reaction rate record
global TIME_M;
global OLD_TIME_M;
global Meta_VEL;
TIME_M=0;
OLD_TIME_M=0;
Meta_VEL=zeros(1,9);

global TIME_N;
global OLD_TIME;
global Gs_VEL;
TIME_N=0;
OLD_TIME=0;
Gs_VEL=zeros(1,9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ini=RAC4leafMetaIni;% Initial values
time=2100;%Simulation time
[Tt, d]=ode15s(@RAC4leafMetaMB,[0,time], Ini); 
global Result;

Result =[Tt,d]; 

vCO2s=Result(:,5).*(Result(:,3)-Result(:,2));
figure;
subplot(2,2,1); plot(Gs_VEL(:,1)/60-5,Gs_VEL(:,2),'k.');title('A');ylim([0,60]);title('A');xlabel('Time (min)');xlim([0,30]);
subplot(2,2,2); plot(Gs_VEL(:,1)/60-5,Gs_VEL(:,9),'k.');title('Transpiration');xlabel('Time (min)');xlim([0,30]);
subplot(2,2,3); plot(Result(:,1)/60-5,Result(:,6),'k.');title('Tleaf');xlabel('Time (min)');xlim([0,30]);
subplot(2,2,4); plot(Result(:,1)/60-5,Result(:,5)*1.6,'k.');title('Gs');xlabel('Time (min)');xlim([0,30]);
if SpeciesCol==1
saveas(gcf,['../Results/','Induction_maize','.tif']);
end
if SpeciesCol==2
saveas(gcf,['../Results/','Induction_sorghum','.tif'])
end
if SpeciesCol==3
saveas(gcf,['../Results/','Induction_sugarcane','.tif'])
end
figure;
plot(Gs_VEL(:,1)/60-5,Gs_VEL(:,11)); title('leakiness');ylim([0,0.6]);xlim([0,30]);


