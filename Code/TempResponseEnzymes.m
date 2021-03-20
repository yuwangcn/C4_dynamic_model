%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TempCorr=TempResponseEnzymes(Temp_leaf)
%WY Temp response 202001
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Tempreature response of enzymes
R=8.3144598;% m2 kg s-2 K-1 mol-1
Ea_KCA=40.9*1000;
dS_KCA=0.21*1000;
Hd_KCA=64.5*51000;
Ea_Vpmax=94.8*1000;
dS_Vpmax=0.25*1000;
Hd_Vpmax=73.3*1000;
Ea_PPDK=58.1*1000;
Ea_Vcmax=78*1000;
Ea_Kc=64.2*1000;
Ea_Ko=10.5*1000;
Vm_OC_25=0.18;
Ea_Vm_OC=55.3*1000;
%Ea_Jmax=43.1*1000;%grow at 25oc Bernacchi 2003
Ea_Jmax=41*1000;%grow at 28oc %% linear correlation between GrowTemp and Ea Bernacchi 2003
%Ea_Jmax=77900;Hd_Jmax=191929;dS_Jmax=627;%MASSAD 2007

Q10=2;

TempCorr_V1=exp(Ea_KCA*(Temp_leaf+273.15-298.15)/(298.15*R*(Temp_leaf+273.15)))*(1+exp((298.15*dS_KCA-Hd_KCA)/(298.15*R)))/(1+exp(((Temp_leaf+273.15)*dS_KCA-Hd_KCA)/((Temp_leaf+273.15)*R)));
TempCorr_V2=exp(Ea_Vpmax *((Temp_leaf+273.15)-298.15)/(298.15*R*(Temp_leaf+273.15)))*(1+exp((298.15*dS_Vpmax-Hd_Vpmax)/(298.15*R)))/(1+exp(((Temp_leaf+273.15)*dS_Vpmax-Hd_Vpmax)/((Temp_leaf+273.15)*R)));
TempCorr_V5=exp(Ea_PPDK*((Temp_leaf+273.15)-298.15)/(298.15*R*(Temp_leaf+273.15)));
TempCorr_V6=exp(Ea_Vcmax*((Temp_leaf+273.15)-298.15)/(298.15*R*(Temp_leaf+273.15)));
TempCorr_KmCO2_6=exp(Ea_Kc*((Temp_leaf+273.15)-298.15)/(298.15*R*(Temp_leaf+273.15)));
TempCorr_KmO2_6=exp(Ea_Ko*((Temp_leaf+273.15)-298.15)/(298.15*R*(Temp_leaf+273.15)));
TempCorr_Vm_OC=exp(Ea_Vm_OC*((Temp_leaf+273.15)-298.15)/(298.15*R*(Temp_leaf+273.15)));
TempCorr_Jmax=exp(Ea_Jmax*((Temp_leaf+273.15)-298.15)/(298.15*R*(Temp_leaf+273.15)));
%TempCorr_Jmax=exp(Ea_Jmax*(Temp_leaf+273.15-298.15)/(298.15*R*(Temp_leaf+273.15)))*(1+exp((298.15*dS_Jmax-Hd_Jmax)/(298.15*R)))/(1+exp(((Temp_leaf+273.15)*dS_Jmax-Hd_Jmax)/((Temp_leaf+273.15)*R)));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sc25=3*10^4;%ubarL/mmol ubarL/mmol
So25=81.58;%kpa L /mmol
%Henry's law constant for CO2 and O2
C_CO2=2400;
C_O2=1700;
kH_CO2=1/Sc25*exp(C_CO2*(1/(Temp_leaf+273.15)-1/298));
kH_O2=1/So25*exp(C_O2*(1/(Temp_leaf+273.15)-1/298));
Sc=1/kH_CO2;
So=1/kH_O2;
TempCorr_KmCO2_6=TempCorr_KmCO2_6*Sc25/Sc;
TempCorr_KmO2_6=TempCorr_KmO2_6*So25/So;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%% ME and MDH and other enzymes
TempCorr_Vm_Enz = Q10^((Temp_leaf-25)/10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TempCorr(1)=TempCorr_V1;
TempCorr(2)=TempCorr_V2;
TempCorr(3)=TempCorr_V5;
TempCorr(4)=TempCorr_V6;
TempCorr(5)=TempCorr_KmCO2_6;
TempCorr(6)=TempCorr_KmO2_6;
TempCorr(7)=TempCorr_Vm_OC;
TempCorr(8)=TempCorr_Jmax;
TempCorr(9)=TempCorr_Vm_Enz;
;
%Chl_O2a=Chl_O2*kH_O2;%WY2018
%Chl_O2=Chl_O2*kH_O2;
% %other enzymes
% Q10_1=1.93;
% Q10_2=2;
% Q10_3=2;
% Q10_5=2;
% Q10_6=2;
% Q10_7=2;
% Q10_8=2;
% Q10_9=2;
% Q10_10=2;
% Q10_13=2;
% Q10_23=2;
% Q10_112=1.81;
% Q10_113=2;
% Q10_121=2;
% Q10_122=2.01;
% Q10_123=2;
% Q10_124=2;
% Q10_131=2;
% Q10_51=2;
% Q10_52=1.60;
% Q10_55=2;
% Q10_56=2;
% Q10_57=2;
% Q10_58=2;
% 
% Ru_Act=-3E-05*Tleaf^3 + 0.0013*Tleaf^2 - 0.0106*Tleaf + 0.8839;%Rubisco activition state
% PsV1 =PsV1_0*Ru_Act*Q10_1^((Tleaf-25)/10);
% PsV2 =PsV2_0*Q10_2^((Tleaf-25)/10);
% PsV3 =PsV3_0*Q10_3^((Tleaf-25)/10);
% PsV5 =PsV5_0*Q10_5^((Tleaf-25)/10);
% PsV6 =PsV6_0*Q10_6^((Tleaf-25)/10);
% PsV7 =PsV7_0*Q10_7^((Tleaf-25)/10);
% PsV8 =PsV8_0*Q10_8^((Tleaf-25)/10);
% PsV9 =PsV9_0*Q10_9^((Tleaf-25)/10);
% PsV10=PsV10_0*Q10_10^((Tleaf-25)/10);
% PsV13=PsV13_0*Q10_13^((Tleaf-25)/10);
% PsV23=PsV23_0*Q10_23^((Tleaf-25)/10);
% PrV112=PrV112_0*Q10_112^((Tleaf-25)/10);
% PrV113=PrV113_0*Q10_113^((Tleaf-25)/10);
% PrV121=PrV121_0*Q10_121^((Tleaf-25)/10);
% PrV122=PrV122_0*Q10_122^((Tleaf-25)/10);
% PrV123=PrV123_0*Q10_123^((Tleaf-25)/10);
% PrV124=PrV124_0*Q10_124^((Tleaf-25)/10);
% PrV131=PrV131_0*Q10_131^((Tleaf-25)/10);
% SUCSV51=SUCSV51_0*Q10_51^((Tleaf-25)/10);
% SUCSV52=SUCSV52_0*Q10_52^((Tleaf-25)/10);
% SUCSV55=SUCSV55_0*Q10_55^((Tleaf-25)/10);
% SUCSV56=SUCSV56_0*Q10_56^((Tleaf-25)/10);
% SUCSV57=SUCSV57_0*Q10_57^((Tleaf-25)/10);
% SUCSV58=SUCSV58_0*Q10_58^((Tleaf-25)/10);
% 
% PrV111= PsV1* 0.24;



