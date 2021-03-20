function LMEnz_v=RAC4LeafMetaVel(t,s)
global gsxsen;
global Tao_MDH;
global Tao_PEPC;
global V6sen;
global RedoxEnyAct;
global EAPPDK;
global kdcon;
global ki;
global kd;
global ainter;
global KValue;
global RatioPPDK;
global Pvpr8M;
%global I;
global TIME_M;
global OLD_TIME_M;
global Meta_VEL;
global Para_mata;
global WeatherTemperature;
global Air_CO2;
global WeatherRH;
global WeatherWind;
global Radiation_PAR;
global Radiation_NIR;
global Radiation_LW;
global OLD_TIME;
global TIME_N;
global Gs_VEL;
global PhiLeaf;
global GsResponse;
global BallBerryInterceptC4;
global BallBerrySlopeC4;
global VmaxC4;
R=8.314472E-3;%Gas constant KJ mole^{-1} K^{-1}
Convert=1E6/(2.35E5); %Convert W m^{-2} to u moles m^{-2} s^{-1}
Boltzman=5.6697E-8; % Stefan-Boltzmann constant W m^{-2} K^{-4}
LatentHeatVaporization=44000.0;%J mole^{-1}
Pressure=101325.0; % Standard atmospheric pressure Pa
ConstantsCp=29.3;
Rd25=0.6;%For steasy state model
PhotosynthesQ10=2;
I=Radiation_PAR*Convert/1000*0.85;
Ci=s(1);
Cb=s(2);
Eb=s(3);
Gs=s(4);
Tleaf=s(5);
H2Oou=s(6);
CO2i=s(7);
LeafTemperature=Tleaf;
Gsw=1.6*Gs;

TempFactor=TempResponseEnzymes(LeafTemperature);
%TempFactor=TempResponseEnzymes(WeatherTemperature);
TempCorr_V1=TempFactor(1);
TempCorr_V2=TempFactor(2);
TempCorr_V5=TempFactor(3);
TempCorr_V6=TempFactor(4);
TempCorr_KmCO2_6=TempFactor(5);
TempCorr_KmO2_6=TempFactor(6);
TempCorr_Vm_OC=TempFactor(7);
TempCorr_Jmax=TempFactor(8);
TempCorr_Vm_Enz=TempFactor(9);

% structure to store the Michaelis-Menten kinetic parameters
KmCO2_1=KValue(1,1);  Ke_1=KValue(1,2);
KmHCO3_2=KValue(2,1);  KmPEP_2=KValue(2,2);   Kimal_2=KValue(2,3);
KmNADPH_3=KValue(3,1);  KmOAA_3=KValue(3,2);  KmNADP_3=KValue(3,3);  Kmmal_3=KValue(3,4);  Ke_3=KValue(3,5);
KmCO2_4=KValue(4,1);  KmNADP_4=KValue(4,2);  KmNADPH_4=KValue(4,3);  KmPyr_4=KValue(4,4);  Kmmal_4=KValue(4,5);  Ke_4=KValue(4,6);
KiPEP_5=KValue(5,1);  KmATP_5=KValue(5,2);  KmPyr_5=KValue(5,3);
KmCO2_6=KValue(6,1)*TempCorr_KmCO2_6;  KmO2_6=KValue(6,2)*TempCorr_KmO2_6;  KmRuBP_6=KValue(6,3);  KiPGA_6=KValue(6,4);  KiFBP_6=KValue(6,5);  KiSBP_6=KValue(6,6);  KiPi_6=KValue(6,7);  KiNADPH_6=KValue(6,8);
KmADP_7=KValue(7,1);  KmATP_7=KValue(7,2);  KmPGA_7=KValue(7,3);
KmDPGA_8=KValue(8,1);  KmNADPH_8=KValue(8,2);
Ke_9=KValue(9,1); 
KmDHAP_10=KValue(10,1);  KmFBP_10=KValue(10,2);  KmGAP_10=KValue(10,3);  Ke_10=KValue(10,4);
KiF6P_11=KValue(11,1);  KiPi_11=KValue(11,2);  KmFBP_11=KValue(11,3);  Ke_11=KValue(11,4);
KmDHAP_12=KValue(12,1);  KmE4P_12=KValue(12,2);  Ke_12=KValue(12,3);
KiPi_13=KValue(13,1);  KmSBP_13=KValue(13,2);  Ke_13=KValue(13,3);
KmE4P_14=KValue(14,1);  KmF6P_14=KValue(14,2);  KmGAP_14=KValue(14,3);  KmXu5P=KValue(14,4);  Ke_14=KValue(14,5);
KmGAP_15=KValue(15,1);  KmRi5P_15=KValue(15,2);  KmS7P_15=KValue(15,3);  KmXu5P_15=KValue(15,4);  Ke_15=KValue(15,5);
Ke_16=KValue(16,1);
Ke_17=KValue(17,1);
KiADP_18=KValue(18,1);  Ki_ADP_18=KValue(18,2);  KiPGA_18=KValue(18,3);  KiPi_18=KValue(18,4);  KiRuBP_18=KValue(18,5);  KmATP_18=KValue(18,6);  KmRu5P_18=KValue(18,7);  Ke_18=KValue(18,8);

KmADP_7Mchl=KValue(19,1);  KmATP_7Mchl=KValue(19,2);  KmPGA_7Mchl=KValue(19,3);
KmDPGA_8Mchl=KValue(20,1);  KmNADPH_8Mchl=KValue(20,1);

KiADP_Starch=KValue(21,1);  KmATP_Starch=KValue(21,2);  KmG1P_Starch=KValue(21,3);  KaF6P_Starch=KValue(21,4);  KaFBP_Starch=KValue(21,5);  KaPGA_Starch=KValue(21,6);  Ke_Starch1=KValue(21,7);   Ke_Starch2=KValue(21,8);
KmPGA_PGASink=KValue(22,1);
KmDHAP_Suc1=KValue(23,1);  KmGAP_Suc1=KValue(23,2);  KmFBP_Suc1=KValue(23,3);  Ke_Suc1=KValue(23,4);
KiF26BP_Suc2=KValue(24,1);  KiF6P_Suc2=KValue(24,2);  KiPi_Suc2=KValue(24,3);  KmFBP_Suc2=KValue(24,4);  Ke_Suc2=KValue(24,5);
Ke_Suc5=KValue(25,1);  Ke_Suc6=KValue(25,2);
KmG1P_Suc7=KValue(26,1);  KmPPi_Suc7=KValue(26,2);  KmUDPG_Suc7=KValue(26,3);  KmUTP_Suc7=KValue(26,4);  Ke_Suc7=KValue(26,5);
KiFBP_Suc8=KValue(27,1);  KiPi_Suc8=KValue(27,2);  KiSuc_Suc8=KValue(27,3);  KiSucP_Suc8=KValue(27,4);  KiUDP_Suc8=KValue(27,5);  KmF6P_Suc8=KValue(27,6);  KmUDPG_Suc8=KValue(27,7);  Ke_Suc8=KValue(27,8); 
KmSuc_Suc9=KValue(28,1);  KmSucP_Suc9=KValue(28,2);  Ke_Suc9=KValue(28,3);
KmSuc_Suc10=KValue(29,1);
KiADP_Suc3=KValue(30,1);  KIDHAP_Suc3=KValue(30,2);  KmATP_Suc3=KValue(30,3);  KmF26BP_Suc3=KValue(30,4);  KmF6P_Suc3=KValue(30,5);  Ke_Suc3=KValue(30,6);
KiF6P_Suc4=KValue(31,1);  KiPi_Suc4=KValue(31,2);  KmF26BP_Suc4=KValue(31,3);
KePi=KValue(36,1);
KmADP_ATPM=KValue(32,1);  KmATP_ATPM=KValue(32,2);  KmPi_ATPM=KValue(32,3);  X=KValue(32,4);  Y=KValue(32,5);  F=KValue(32,6);  Q=KValue(32,7);  D=KValue(32,8); Ke_ATPM=KValue(32,9);
KmNADP_NADPHM=KValue(33,1);  KmNADPH_NADPHM=KValue(33,2); Ke_NADPHM=KValue(33,3);  E=KValue(33,4);
KmADP_ATPB=KValue(34,1);  KmPi_ATPB=KValue(34,2);   KmATP_ATPB=KValue(34,3);  Ke_ATPB=KValue(34,4);   G=KValue(34,5);
KmNADP_NADPHB=KValue(37,1);  KmNADPH_NADPHB=KValue(37,2); Ke_NADPHB=KValue(37,3);
Voaa=KValue(35,1);  Vmal=KValue(35,2);  Vpyr=KValue(35,3);  Vpep=KValue(35,4);  Vt=KValue(35,5);  Vleak=KValue(35,6); Vpga=KValue(35,7);
KmCO2_PR1=KValue(38,1)*TempCorr_KmCO2_6; KmO2_PR1=KValue(38,2)*TempCorr_KmO2_6;  KmRuBP_PR1=KValue(38,3);  KiPGA_PR1=KValue(38,4);  KiFBP_PR1=KValue(38,5);  KiSBP_PR1=KValue(38,6);  KiPi_PR1=KValue(38,7);  KiNADPH_PR1=KValue(38,8);
KmPGCA_PR2=KValue(39,1);  KiPI_PR2=KValue(39,2);  KiGCA_PR2=KValue(39,3);
KmGCA_PR3=KValue(40,1);
Ke_PS4=KValue(41,1);  KmGOA_PS4=KValue(41,2);  KmGLU_PS4=KValue(41,3);  KiGLY_PS4=KValue(41,4);
KmGLY_PS5=KValue(42,1);  KiSER_PS5=KValue(42,2);
Ke_PR6=KValue(43,1);  KmGOA_PR6=KValue(43,2);  KmSER_PR6=KValue(43,3);  KmGLY_PR6=KValue(43,4);
Ke_PR7=KValue(44,1);  KiHPR_PR7=KValue(44,2);  KmHPR_PR7=KValue(44,3);
Ke_PR8=KValue(45,1);  KmATP_PR8=KValue(45,2);  KmGCEA_PR8=KValue(45,3);  KiPGA_PR8=KValue(45,4);
KmGCA_PR9=KValue(46,1);  KiGCEA_PR9=KValue(46,2);
KmGCEA_PR10=KValue(47,1);  KiGCA_PR10=KValue(47,2);
KmPGA_62=KValue(48,1); KmPEP_62=KValue(48,2); Ke_62=KValue(48,3);


MC_HCO3= s(1+7);
MC_OAA=s(2+7);
MC_PEP=s(3+7);
MC_malate=s(4+7);
MC_pyruvate=s(5+7);
MC_PGA=s(6+7);
MC_FBP=s(7+7);
MC_UDPG=s(8+7);
MC_SUCP=s(9+7);
MC_SUC=s(10+7);
MC_F26BP=s(11+7);
MC_ATP=s(12+7);
MC_T3P=s(13+7);
MC_HexP=s(14+7);
MC_Sucrose=s(15+7);
Mchl_OAA= s(16+7);
Mchl_malate =s(17+7);
Mchl_PEP =s(18+7);
Mchl_pyruvate= s(19+7);
Mchl_NADPH= s(20+7);
Mchl_ATP= s(21+7);
Mchl_PGA= s(22+7);
Mchl_DPGA= s(23+7);
Mchl_T3P= s(24+7);
BSC_T3P= s(25+7);
BSC_PGA= s(26+7);
BSC_malate= s(27+7);
BSC_pyruvate= s(28+7);
BSC_CO2=s(29+7);
Bchl_CO2= s(30+7);
Bchl_RuBP= s(31+7);
Bchl_PGA= s(32+7);
Bchl_DPGA= s(33+7);
Bchl_ATP=s(34+7);
Bchl_NADPH= s(35+7);
Bchl_SBP= s(36+7);
Bchl_S7P= s(37+7);
Bchl_FBP= s(38+7);
Bchl_E4P= s(39+7);
Bchl_Starch= s(40+7);
Bchl_Rubisco= s(41+7);
Bchl_T3P= s(42+7);
Bchl_HexP= s(43+7);
Bchl_Pent =s(44+7);
Bchl_malate= s(45+7);
Bchl_pyruvate= s(46+7);
Bchl_PGCA=s(47+7);
Bchl_GCA=s(48+7);
Bchl_GCEA=s(49+7);
Bper_GCA=s(50+7);
Bper_GOA=s(51+7);
Bper_GLY=s(52+7);
Bper_SER=s(53+7);
Bper_HPR=s(54+7);
Bper_GCEA=s(55+7);
MC_CO2=s(56+7);
Bchl_PPi=s(57+7);
Bchl_ADPG=s(58+7);
MC_Glu=s(59+7);
MC_OxoG=s(60+7);
MC_Asp=s(61+7);
MC_Ala=s(62+7);
BSC_OxoG=s(63+7);
BSC_Glu=s(64+7);
BSC_Asp=s(65+7);
BSC_Ala=s(66+7);
BSC_OAA=s(67+7);
BSC_PEP=s(68+7);
BSC_ATP=s(69+7);
Bchl_OAA=s(70+7);

MC_O2=s(71+7);
Mchl_O2=s(72+7);
BSC_O2=s(73+7);
Bchl_O2=s(74+7);
Bchl_PEP=s(75+7);%%%%%%%WY PPDK in BSCytosol
Mchl_GCEA=s(76+7);
Bmito_OAA=s(77+7);
Bmito_MAL=s(78+7);
Bmito_PYR=s(79+7);
Bmito_CO2=s(80+7);
Bmito_NADH=s(81+7);
Bchl_Asp=s(82+7);
Bchl_Ala=s(83+7);
Mchl_Asp=s(84+7);
Mchl_Ala=s(85+7);
E_PPDK_Mchl=s(86+7);% E_PPDK active
EP_PPDK_Mchl=s(87+7);%EP_PPDK inactive

global Velocity_s;
 
Vm_1=Velocity_s(1)*TempCorr_V1;
Vm_2=Velocity_s(2)*TempCorr_V2;          
Vm_3=Velocity_s(3)*TempCorr_Vm_Enz;
Vm_4=Velocity_s(4)*TempCorr_Vm_Enz;        
Vm_5=Velocity_s(5)*TempCorr_V5;
Vm_6=Velocity_s(6)*TempCorr_V6;
Vm_78=Velocity_s(7)*TempCorr_Vm_Enz;
Vm_10=Velocity_s(9)*TempCorr_Vm_Enz;
Vm_11=Velocity_s(10)*TempCorr_Vm_Enz;
Vm_12=Velocity_s(11)*TempCorr_Vm_Enz;
Vm_13=Velocity_s(12)*TempCorr_Vm_Enz;
Vm_14=Velocity_s(13)*TempCorr_Vm_Enz;
Vm_15=Velocity_s(14)*TempCorr_Vm_Enz;
Vm_18=Velocity_s(15)*TempCorr_Vm_Enz;
Vm_78Mchl=Velocity_s(16)*TempCorr_Vm_Enz;
Vm_Starch=Velocity_s(18)*TempCorr_Vm_Enz;
Vm_PGASink=Velocity_s(19)*TempCorr_Vm_Enz;
Vm_Suc1=Velocity_s(20)*TempCorr_Vm_Enz;
Vm_Suc2=Velocity_s(21)*TempCorr_Vm_Enz;
Vm_Suc7=Velocity_s(22)*TempCorr_Vm_Enz;
Vm_Suc8=Velocity_s(23)*TempCorr_Vm_Enz;
Vm_Suc9=Velocity_s(24)*TempCorr_Vm_Enz;
Vm_Suc10=Velocity_s(25)*TempCorr_Vm_Enz;
Vm_Suc3=Velocity_s(26)*TempCorr_Vm_Enz;
Vm_Suc4=Velocity_s(27)*TempCorr_Vm_Enz;
Jmax=Velocity_s(29)*TempCorr_Jmax;
Vm_ATPM=Velocity_s(30)*TempCorr_Vm_Enz;
Vm_NADPHM=Velocity_s(31)*TempCorr_Vm_Enz; 
Vm_ATPB=Velocity_s(32)*TempCorr_Vm_Enz;
Vm_NADPHB=Velocity_s(33)*TempCorr_Vm_Enz;
%Vm_PR1=Vm_6*0.11*TempCorr_Vm_OC;
Vm_PR2=Velocity_s(35)*TempCorr_Vm_Enz;
Vm_PR3=Velocity_s(36)*TempCorr_Vm_Enz;
Vm_PR4=Velocity_s(37)*TempCorr_Vm_Enz;
Vm_PR5=Velocity_s(38)*TempCorr_Vm_Enz;
Vm_PR6=Velocity_s(39)*TempCorr_Vm_Enz;
Vm_PR7=Velocity_s(40)*TempCorr_Vm_Enz;
Vm_PR8=Velocity_s(41)*TempCorr_Vm_Enz;
VTgca_PR9=Velocity_s(42)*TempCorr_Vm_Enz;
VTgcea_PR10=Velocity_s(43)*TempCorr_Vm_Enz;
Vm_62=Velocity_s(44)*TempCorr_Vm_Enz;
Vtp_Bchl=Velocity_s(45)*TempCorr_Vm_Enz;
Vtp_Mchl=Velocity_s(46)*TempCorr_Vm_Enz;
Vm_Sta1=Velocity_s(47)*TempCorr_Vm_Enz;
Vm_Sta2=Velocity_s(48)*TempCorr_Vm_Enz;
Vm_Sta3=Velocity_s(49)*TempCorr_Vm_Enz;
Vm_OAA_M=Velocity_s(50)*TempCorr_Vm_Enz;
Vm_PYR_B=Velocity_s(51)*TempCorr_Vm_Enz;
Vm_PYR_M=Velocity_s(52)*TempCorr_Vm_Enz;
Vm_PEP_M=Velocity_s(53)*TempCorr_Vm_Enz;
Pmal=Velocity_s(54)*TempCorr_Vm_Enz;
Ppyr=Velocity_s(55)*TempCorr_Vm_Enz;
Pco2=Velocity_s(56)*TempCorr_Vm_Enz;
PC3P=Velocity_s(57)*TempCorr_Vm_Enz;
Pco2_B=Velocity_s(58)*TempCorr_Vm_Enz;

Vm_MAL_B=Vm_PEP_M*TempCorr_Vm_Enz;
Vm_MAL_M=Vm_PEP_M*TempCorr_Vm_Enz;
gm=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Bchl_CA;%assume constant concentration
global Bchl_CN;
global Bchl_CP;
global MC_CU;
global MC_CA;
global MC_CP;
global MC_UTP;
global Mchl_CP;
global Mchl_CA;
global Mchl_CN;
Mchl_NADP = Mchl_CN-Mchl_NADPH;
Mchl_Pi = Mchl_CP-Mchl_PGA-2*Mchl_DPGA-Mchl_T3P-Mchl_ATP-Mchl_PEP;
Mchl_GAP = Ke_9*Mchl_T3P/(1+Ke_9);
Mchl_DHAP = Mchl_T3P/(1+Ke_9);
Mchl_ADP = Mchl_CA-Mchl_ATP;

MC_UDP = MC_CU-MC_UTP-MC_UDPG;
MC_PiT = MC_CP-2*MC_FBP-2*MC_F26BP-MC_PGA-MC_T3P-MC_HexP-MC_SUCP-MC_UTP-MC_ATP-MC_PEP;%1
MC_Pi = (sqrt(KePi^2+4*KePi*MC_PiT)-KePi)/2;
MC_PPi = MC_PiT-MC_Pi;
MC_GAP = Ke_9*MC_T3P/(1+Ke_9);
MC_DHAP = MC_T3P/(1+Ke_9);
MC_G6P = MC_HexP/(1/Ke_Suc5+Ke_Suc6+1);
MC_G1P = Ke_Suc6*MC_HexP/(1/Ke_Suc5+Ke_Suc6+1);
MC_F6P = (MC_HexP/Ke_Suc5)/(1/Ke_Suc5+Ke_Suc6+1);
MC_ADP = MC_CA-MC_ATP;
BSC_GAP = Ke_9*BSC_T3P/(1+Ke_9);
BSC_DHAP = BSC_T3P/(1+Ke_9);
Bchl_NADP = Bchl_CN-Bchl_NADPH;
Bchl_Xu5P = (Bchl_Pent/Ke_17)/(1/Ke_16+1/Ke_17+1);
Bchl_Ru5P = Bchl_Pent/(1/Ke_16+1/Ke_17+1);
Bchl_Ri5P = (Bchl_Pent/Ke_16)/(1/Ke_16+1/Ke_17+1);

Bchl_Pi = Bchl_CP-Bchl_PGA-2*Bchl_DPGA-Bchl_T3P-2*Bchl_FBP-Bchl_HexP-Bchl_E4P-2*Bchl_SBP-Bchl_S7P-Bchl_Pent-2*Bchl_RuBP-Bchl_ATP-Bchl_PGCA-Bchl_PPi-Bchl_PEP;%%%WY
global Pim;
Pim=Bchl_Pi;
Bchl_GAP = Ke_9*Bchl_T3P/(1+Ke_9);
Bchl_DHAP = Bchl_T3P/(1+Ke_9);
Bchl_G6P = Bchl_HexP/(1/Ke_Starch1+Ke_Starch2+1);
Bchl_G1P = Ke_Starch2*Bchl_HexP/(1/Ke_Starch1+Ke_Starch2+1);
Bchl_F6P = (Bchl_HexP/Ke_Starch1)/(1/Ke_Starch1+Ke_Starch2+1);
Bchl_ADP = Bchl_CA-Bchl_ATP-Bchl_ADPG;
Bmito_NAD=1-Bmito_NADH;

global CI;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Enzyme activation %WY 201902
%PPDK
global PPDKRP;
AMP=0;
PPi=0;
ET_PKRP=PPDKRP;%0.000025*3;%%%change 0.0001
Kcat_EA_PPDKRP_I= 1.125; % s-1
Kcat_EA_PPDKRP_A= 0.578; % s-1
% inactivation
Km_EA_PPDKRP_I_ADP=0.052; %mM
Ki_EA_PPDKRP_I_Pyr=0.08;
Km_EA_PPDKRP_I_E=0.0012;
% activation
Km_EA_PPDKRP_A_Pi= 0.65;
Ki_EA_PPDKRP_A_AMP=0.4;
Km_EA_PPDKRP_A_EP=0.0007;
Ki_EA_PPDKRP_A_ADP=0.085;
Ki_EA_PPDKRP_A_PPI=0.16;
% inactivation rate
vEA_PPDKRP_I =ET_PKRP*Kcat_EA_PPDKRP_I* E_PPDK_Mchl*Mchl_ADP/ ((E_PPDK_Mchl+Km_EA_PPDKRP_I_E)*(Mchl_ADP+ Km_EA_PPDKRP_I_ADP*(1+Mchl_pyruvate/ Ki_EA_PPDKRP_I_Pyr)));
% activation rate
vEA_PPDKRP_A= ET_PKRP*Kcat_EA_PPDKRP_A* EP_PPDK_Mchl*Mchl_Pi/(( EP_PPDK_Mchl+ Km_EA_PPDKRP_A_EP*(1+Mchl_ADP/ Ki_EA_PPDKRP_A_ADP+PPi/ Ki_EA_PPDKRP_A_PPI))*( Mchl_Pi+ Km_EA_PPDKRP_A_Pi*(1+AMP/ Ki_EA_PPDKRP_A_AMP)));
A_PPDK=E_PPDK_Mchl/(E_PPDK_Mchl+EP_PPDK_Mchl);
if EAPPDK==1
Vm_5=Vm_5*A_PPDK;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mchl_ActATPsynthase=s(99);
Mchl_ActGAPDH=s(100);
Mchl_ActNADPMDH=s(101);
Bchl_ActATPsynthase=s(102);
Mchl_ActPEPC=s(103);
Bchl_ActGAPDH=s(104);
Bchl_ActFBPase=s(105);
Bchl_ActSBPase=s(106);
Bchl_ActPRK=s(107);
Bchl_ActRubisco=s(108);
Bchl_ActRCA=s(109);
%enzyme activation 
AACCC=0.05;
ActPEPC0=0.0017*I*1000/0.85+0.05;
ActFBPase0=0.0017*I*1000/0.85+AACCC;
ActSBPase0=0.0017*I*1000/0.85+AACCC;
ActGAPDH0=0.0017*I*1000/0.85+AACCC;
ActPRK0=0.0017*I*1000/0.85/+AACCC;
ActATPsyn0=0.0017*I*1000/0.85+AACCC;
ActNADPMDH0=0.0017*I*1000/0.85+0.05;

if ActPEPC0>1
  ActPEPC0=1;
end
if ActFBPase0>1
  ActFBPase0=1;
end
if ActSBPase0>1
  ActSBPase0=1;
end
if ActGAPDH0>1
  ActGAPDH0=1;
end
if ActPRK0>1
  ActPRK0=1;
end
if ActATPsyn0>1
  ActATPsyn0=1;
end
if ActNADPMDH0>1
  ActNADPMDH0=1;
end

global taoRub;
tao_ActPEPC =60*Tao_PEPC;
tao_ActFBPase =1.878*60;%1.878*60;
tao_ActSBPase =3.963*60;%3.963*60;
tao_ActATPsynthase=0.5*60;
tao_ActGAPDH=1*60/10;
tao_ActPRK=1*60/10;
tao_ActNADPMDH=0.965*60*Tao_MDH;%0.5*60*4*Tao_MDH;
KaRac=12.4;%mg m-2
tao_ActRubisco=taoRub*60;%2.5;%KTaoRac/Rac; % s

% kcat_ATPsyn=16;
% kcat_NADPMDH=1520;
% kcat_PEPC=66;
% kcat_Rubisco=4.1;
% Kcat_GAPDH=50;
% kcat_FBP=22.9;
% kcat_SBP=81;
% kcat_PRK=615;

vATPsynthase_Act_Mchl=(ActATPsyn0-Mchl_ActATPsynthase)*1/tao_ActATPsynthase;%*Vm_ATPM/kcat_ATPsyn;
vNADPMDH_Act=(ActNADPMDH0-Mchl_ActNADPMDH)*1/tao_ActNADPMDH;%*Vm_3/kcat_NADPMDH;
vGAPDH_Act_Mchl=(ActGAPDH0-Mchl_ActGAPDH)*1/tao_ActGAPDH;
vATPsynthase_Act_Bchl=(ActATPsyn0-Bchl_ActATPsynthase)*1/tao_ActATPsynthase;%;
vPEPC_Act=(ActPEPC0-Mchl_ActPEPC)*1/tao_ActPEPC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%WY 202010
ActRca0=0.0017*I*1000/0.85+0.06;
if ActRca0>1
  ActRca0=1;
end
tao_ActRca =0.7594*60/10;
vRCA_Act=(ActRca0-Bchl_ActRCA)*1/tao_ActRca;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vGAPDH_Act_Bchl=(ActGAPDH0-Bchl_ActGAPDH)*1/tao_ActGAPDH;
vFBPase_Act=(ActFBPase0-Bchl_ActFBPase)*1/tao_ActFBPase;
vSBPase_Act=(ActSBPase0-Bchl_ActSBPase)*1/tao_ActSBPase;
vPRK_Act=(ActPRK0-Bchl_ActPRK)*1/tao_ActPRK;


Vm_ATPM=Mchl_ActATPsynthase*Vm_ATPM;

if RedoxEnyAct==1
Vm_3=Mchl_ActNADPMDH*Vm_3;
Vm_2=Mchl_ActPEPC*Vm_2;
Vm_ATPB=Bchl_ActATPsynthase*Vm_ATPB;
Vm_78Mchl=(Mchl_ActGAPDH)*Vm_78Mchl;
Vm_78=(Bchl_ActGAPDH)*Vm_78;
Vm_11=(Bchl_ActFBPase)*Vm_11;
Vm_13=(Bchl_ActSBPase)*Vm_13;
Vm_18=(Bchl_ActPRK)*Vm_18;
end
global PRac
Rac=216.9/tao_ActRubisco*60;%%216.9 min mg m-2
ActRubisco0=Rac*Bchl_ActRCA/(KaRac+Rac*Bchl_ActRCA);
vRubisco_Act=(ActRubisco0-Bchl_ActRubisco)*1/tao_ActRubisco;%;
if PRac==1
Vm_6=(Bchl_ActRubisco)*Vm_6*V6sen;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Bper_GLU;
global Bper_KG;
global Bper_NADH;
global Bper_NAD;

Sc=3*10^4;%3.36*10^4;%ubarL/mmolubarL/mmol
vs=20;% m2/L

if Para_mata==0
vinf=gm*Sc*10^(-3)*(CI-MC_CO2);
end
if Para_mata==1
vinf=gm*Sc*10^(-3)*(Ci/(3 * 10^4)-MC_CO2);
end

%C4 Cycle 5
global v1;
v1=Vm_1*(MC_CO2-MC_HCO3/Ke_1)/(KmCO2_1+MC_CO2);

global v2;
Kimal_2=15;
Kimal_2n=1.5;
Vm_2P=Mchl_ActPEPC*Vm_2;
Vm_2OH=Vm_2-Vm_2P;
v_2P=Vm_2*MC_HCO3*MC_PEP/(MC_PEP+KmPEP_2*(1+MC_malate/Kimal_2+MC_OAA/1))/(MC_HCO3+KmHCO3_2);
v_2OH=Vm_2OH*MC_HCO3*MC_PEP/(MC_PEP+KmPEP_2*(1+MC_malate/(Kimal_2n)+MC_OAA/1))/(MC_HCO3+KmHCO3_2);
v2=v_2P+v_2OH;

global v3;
v3=Vm_3*(Mchl_OAA*Mchl_NADPH-Mchl_NADP*Mchl_malate/Ke_3)/( KmOAA_3* KmNADPH_3*(1+Mchl_OAA/KmOAA_3+ Mchl_NADPH/KmNADPH_3+ Mchl_NADP/KmNADP_3+ Mchl_malate/Kmmal_3+ Mchl_OAA*Mchl_NADPH/(KmOAA_3* KmNADPH_3)+ Mchl_NADP*Mchl_malate/(KmNADP_3* Kmmal_3)));

global v4;
v4=Vm_4*(Bchl_malate*Bchl_NADP-Bchl_pyruvate*Bchl_NADPH*Bchl_CO2/Ke_4)/(Kmmal_4*KmNADP_4)/(1+Bchl_malate/Kmmal_4+Bchl_NADP/KmNADP_4+Bchl_pyruvate/KmPyr_4+Bchl_NADPH/KmNADPH_4+Bchl_CO2/KmCO2_4+Bchl_malate*Bchl_NADP/(Kmmal_4*KmNADP_4)+Bchl_pyruvate*Bchl_NADPH/(KmPyr_4*KmNADPH_4)+Bchl_pyruvate*Bchl_CO2/(KmPyr_4*KmCO2_4)+Bchl_NADPH*Bchl_CO2/(KmNADPH_4*KmCO2_4)+Bchl_pyruvate*Bchl_NADPH*Bchl_CO2/(KmPyr_4*KmNADPH_4*KmCO2_4));
global v5;
v5=Vm_5*Mchl_pyruvate*Mchl_ATP/(Mchl_pyruvate+KmPyr_5*(1+Mchl_PEP/KiPEP_5))/(Mchl_ATP+KmATP_5);
%%%% PPDK in BSC
global v5B;
Vm_5B=RatioPPDK*Vm_5;
v5B=Vm_5B*Bchl_pyruvate*Bchl_ATP/(Bchl_pyruvate+KmPyr_5*(1+Bchl_PEP/KiPEP_5))/(Bchl_ATP+KmATP_5);
%%%%%
% % % Calvin cycle 
global v6;
v6=Vm_6*Bchl_RuBP*Bchl_CO2/((Bchl_CO2+KmCO2_6*(1+Bchl_O2/KmO2_6))*(Bchl_RuBP+KmRuBP_6*(1+Bchl_PGA/KiPGA_6+Bchl_FBP/KiFBP_6+Bchl_SBP/KiSBP_6+Bchl_Pi/ KiPi_6)));
v7=0;%Not used 
v8=0;%Not used 

global v78;
KmPGA_78=5;
KmATP_78=0.3;
KmNADPH_78=0.1;
KmADP_78=0.5;
KmNADP_78=0.5;
%v78=Vm_78*(Bchl_PGA*Bchl_ATP*Bchl_NADPH)/((Bchl_PGA+KmPGA_78)*(Bchl_ATP+KmATP_78*(1+Bchl_ADP/KmADP_78))*(Bchl_NADPH+KmNADPH_78*(1+Bchl_NADP/KmNADP_78)));
v78=Vm_78*(Bchl_PGA*Bchl_ATP*Bchl_NADPH)/((Bchl_PGA+KmPGA_78)*(Bchl_ATP+KmATP_78)*(Bchl_NADPH+KmNADPH_78));
global v10;
global v11;
global v12;
global v13;
global v14;
global v15;
global v18;
v10=Vm_10*(Bchl_GAP*Bchl_DHAP-Bchl_FBP/Ke_10)/(KmGAP_10*KmDHAP_10*(1+Bchl_GAP/KmGAP_10+Bchl_DHAP/KmDHAP_10+Bchl_FBP/KmFBP_10+Bchl_GAP*Bchl_DHAP/(KmGAP_10*KmDHAP_10)));
v11=Vm_11*(Bchl_FBP-Bchl_F6P*Bchl_Pi/Ke_11)/(Bchl_FBP+KmFBP_11*(1+Bchl_F6P/KiF6P_11+Bchl_Pi/KiPi_11));
v12=Vm_12*(Bchl_DHAP*Bchl_E4P-Bchl_SBP/Ke_12)/((Bchl_E4P+KmE4P_12)*(Bchl_DHAP+KmDHAP_12));
v13=Vm_13*(Bchl_SBP-Bchl_Pi*Bchl_S7P/Ke_13)/(Bchl_SBP+KmSBP_13*(1+Bchl_Pi/KiPi_13));
v14=Vm_14*(Bchl_F6P*Bchl_GAP-Bchl_Xu5P*Bchl_E4P/Ke_14)/((Bchl_F6P + KmF6P_14*(1+Bchl_Xu5P/KmXu5P+Bchl_E4P/KmE4P_14))*(Bchl_GAP+KmGAP_14));
v15=Vm_15*(Bchl_GAP*Bchl_S7P-Bchl_Ri5P*Bchl_Xu5P/Ke_15)/((Bchl_GAP+KmGAP_15*(1+Bchl_Xu5P/KmXu5P_15+Bchl_Ri5P/KmRi5P_15))*(Bchl_S7P+KmS7P_15));
v18=Vm_18*(Bchl_ATP*Bchl_Ru5P-Bchl_ADP*Bchl_RuBP/Ke_18)/((Bchl_ATP*(1+Bchl_ADP/KiADP_18)+KmATP_18*(1+Bchl_ADP/Ki_ADP_18))*(Bchl_Ru5P+KmRu5P_18*(1+Bchl_PGA/KiPGA_18+Bchl_RuBP/KiRuBP_18+Bchl_Pi/KiPi_18)));
v7Mchl=0;%Not used;
v8Mchl=0;%Not used ;
global v78Mchl;
KmPGA_78Mchl=5;
KmATP_78Mchl=0.3;
KmNADPH_78Mchl=0.1;
KmADP_78Mchl=0.5;
KmNADP_78Mchl=0.5;

v78Mchl=Vm_78Mchl*(Mchl_PGA*Mchl_ATP*Mchl_NADPH)/((Mchl_PGA+KmPGA_78Mchl)*(Mchl_ATP+KmATP_78Mchl)*(Mchl_NADPH+KmNADPH_78Mchl));
vStarch1=Vm_Starch*Bchl_G1P*Bchl_ATP/((Bchl_G1P+KmG1P_Starch)*((1+Bchl_ADP/KiADP_Starch)*(Bchl_ATP+KmATP_Starch)+KmATP_Starch*Bchl_Pi/(KaPGA_Starch*Bchl_PGA+KaF6P_Starch*Bchl_F6P+KaFBP_Starch*Bchl_FBP)));
% KmATP_Starch2=1;
% Vm_Starch2=0.1/20;
% KmStarch_Starch2=1;
vStarch2=0;%Not used 

%  2.7.7.27
KaPGA_Sta1=0.2;%0.2252;
KmG1P_Sta1=0.06;%0.06;
KmATP_Sta1=0.12;%0.18
KIAPi_ATP_Sta1=2.96;
KmPPi_Sta1=0.033;
KICPP1_ATP_Sta1=13.8E-4;
KmADPG_Sta1=0.24;
KIAADP_ATP_Sta1=2.0;
Ke_Sta1=1.1;

global vSta1;
Vm_Sta1=Vm_Sta1*Bchl_PGA/(Bchl_PGA+KaPGA_Sta1);
vSta1=Vm_Sta1*(Bchl_G1P*Bchl_ATP-Bchl_ADPG*Bchl_PPi/Ke_Sta1)/(KmG1P_Sta1*KmATP_Sta1*(1+Bchl_ADP/KIAADP_ATP_Sta1+Bchl_PPi/KICPP1_ATP_Sta1+Bchl_Pi/KIAPi_ATP_Sta1)*(1+Bchl_G1P/KmG1P_Sta1+Bchl_ATP*(1+Bchl_Pi/KIAPi_ATP_Sta1+Bchl_ADP/KIAADP_ATP_Sta1)/(KmATP_Sta1*(1+Bchl_ADP/KIAADP_ATP_Sta1+Bchl_PPi/KICPP1_ATP_Sta1+Bchl_Pi/KIAPi_ATP_Sta1))+Bchl_ADPG/KmADPG_Sta1+Bchl_PPi/KmPPi_Sta1+Bchl_G1P*Bchl_ATP*(1+Bchl_Pi/KIAPi_ATP_Sta1+Bchl_ADP/KIAADP_ATP_Sta1)/(KmG1P_Sta1*KmATP_Sta1*(1+Bchl_ADP/KIAADP_ATP_Sta1+Bchl_PPi/KICPP1_ATP_Sta1+Bchl_Pi/KIAPi_ATP_Sta1))+Bchl_ADPG*Bchl_PPi/(KmADPG_Sta1*KmPPi_Sta1)));
%  3.6.1.1
KmPPi_Sta2=0.154;
Ke_Sta2=15700.0;
vSta2=Vm_Sta2*(Bchl_PPi-Bchl_Pi*Bchl_Pi/Ke_Sta2)/(Bchl_PPi+KmPPi_Sta2);

%  2.4.1.21
KmADPG_Sta3=0.077;
vSta3=Vm_Sta3*Bchl_ADPG/(Bchl_ADPG+KmADPG_Sta3);
global vStarch;
vStarch=vSta3;


Kmpi_hexp=1.5;
Kmhexp_hexp=1;%1/5;%WY20/02
%KmATP_hexp=0.8;
Vm_Hep=0.0005;
if Bchl_HexP>=0&& Bchl_HexP <6
Vm_Hep=(6-Bchl_HexP)/6*0.0005;
end
if Bchl_HexP>6
 Vm_Hep=0;
end

vhexp=Vm_Hep*(Bchl_Pi/(Kmpi_hexp*(1+Bchl_HexP/Kmhexp_hexp)+Bchl_Pi));
if Bchl_HexP>6
    vhexp=0;
end
global vHexP;
vHexP=vhexp;
vPGASink=Vm_PGASink*MC_PGA/(MC_PGA+KmPGA_PGASink);%MC
global vpgasink;
vpgasink=vPGASink;
global vSuc1;
vSuc1=Vm_Suc1*(MC_GAP*MC_DHAP-MC_FBP/Ke_Suc1)/(KmGAP_Suc1*KmDHAP_Suc1*(1+MC_GAP/KmGAP_Suc1+MC_DHAP/KmDHAP_Suc1+MC_FBP/KmFBP_Suc1+MC_GAP*MC_DHAP/(KmGAP_Suc1*KmDHAP_Suc1)));
vSuc2=Vm_Suc2*(MC_FBP-MC_F6P*MC_Pi/Ke_Suc2)/(KmFBP_Suc2*(1+MC_F26BP/KiF26BP_Suc2)*(1+MC_FBP/(KmFBP_Suc2*(1+MC_F26BP/KiF26BP_Suc2))+MC_Pi/KiPi_Suc2+MC_F6P/KiF6P_Suc2+MC_Pi*MC_F6P/(KiPi_Suc2*KiF6P_Suc2)));
vSuc7=Vm_Suc7*(MC_UTP*MC_G1P-MC_UDPG*MC_PPi/Ke_Suc7)/(KmUTP_Suc7*KmG1P_Suc7*(1+MC_UTP/KmUTP_Suc7+MC_G1P/KmG1P_Suc7+MC_UDPG/KmUDPG_Suc7+MC_PPi/KmPPi_Suc7+MC_UTP*MC_G1P/(KmUTP_Suc7*KmG1P_Suc7)+MC_UDPG*MC_PPi/(KmUDPG_Suc7*KmPPi_Suc7)));
vSuc8=Vm_Suc8*(MC_F6P*MC_UDPG-MC_SUCP*MC_UDP/Ke_Suc8)/((MC_F6P+KmF6P_Suc8*(1+MC_FBP/KiFBP_Suc8))*(MC_UDPG+KmUDPG_Suc8*(1+MC_UDP/KiUDP_Suc8)*(1+MC_SUCP/KiSucP_Suc8)*(1+MC_SUC/KiSuc_Suc8)*(1+MC_Pi/KiPi_Suc8)));
vSuc9=Vm_Suc9*(MC_SUCP-MC_SUC*MC_Pi/Ke_Suc9)/(MC_SUCP+KmSucP_Suc9*(1+MC_SUC/KmSuc_Suc9));
vSuc10=Vm_Suc10*MC_SUC/(MC_SUC+KmSuc_Suc10);
global vSuc;
vSuc= vSuc10;
vSuc3=Vm_Suc3*(MC_ATP*MC_F6P-MC_ADP*MC_F26BP/Ke_Suc3)/((MC_F6P+KmF6P_Suc3*(1+MC_F26BP/KmF26BP_Suc3)*(1+MC_DHAP/KIDHAP_Suc3))*(MC_ATP+KmATP_Suc3*(1+MC_ADP/KiADP_Suc3)));
vSuc4=Vm_Suc4*MC_F26BP/(KmF26BP_Suc4*(1+MC_F26BP/KmF26BP_Suc4)*(1+MC_Pi/KiPi_Suc4)*(1+MC_F6P/KiF6P_Suc4));

%ATP&NADPH 3
global ETRa;
global ETRn;
ETRa=D*((F/2*X*I+Y*Jmax-sqrt((F/2*X*I+Y*Jmax)^2-4*Q*F/2*X*I*Y*Jmax))/(2*Q));
ETRn=E*((F/2*X*I+Y*Jmax-sqrt((F/2*X*I+Y*Jmax)^2-4*Q*F/2*X*I*Y*Jmax))/(2*Q));
vATPM=min(Vm_ATPM,ETRa)*(Mchl_ADP*Mchl_Pi-Mchl_ATP/(Ke_ATPM))/(KmADP_ATPM*KmPi_ATPM*(1+Mchl_ADP/KmADP_ATPM+Mchl_Pi/KmPi_ATPM+Mchl_ATP/KmATP_ATPM+Mchl_ADP*Mchl_Pi/(KmADP_ATPM*KmPi_ATPM)));
vNADPHM=min(Vm_NADPHM,ETRn)*(Mchl_NADP-Mchl_NADPH/Ke_NADPHM)/(KmNADP_NADPHM*(1+Mchl_NADP/KmNADP_NADPHM+Mchl_NADPH/KmNADPH_NADPHM));
% %WY202003
% vATPM=min(Vm_ATPM*(Mchl_ADP*Mchl_Pi-Mchl_ATP/(Ke_ATPM))/(KmADP_ATPM*KmPi_ATPM*(1+Mchl_ADP/KmADP_ATPM+Mchl_Pi/KmPi_ATPM+Mchl_ATP/KmATP_ATPM+Mchl_ADP*Mchl_Pi/(KmADP_ATPM*KmPi_ATPM))),ETRa);
% vNADPHM=min(Vm_NADPHM*(Mchl_NADP-Mchl_NADPH/Ke_NADPHM)/(KmNADP_NADPHM*(1+Mchl_NADP/KmNADP_NADPHM+Mchl_NADPH/KmNADPH_NADPHM)),ETRn);


global vO2_Mchl;
vO2_Mchl=vNADPHM/2;

global U;
global V;
global ETRab
global ETRabl;
global ETRnbl;

ETRab=G*((F*(1-U)*(1-X)*I+(1-V)*(1-Y)*Jmax-sqrt((F*(1-U)*(1-X)*I+(1-V)*(1-Y)*Jmax)^2-4*Q*F*(1-U)*(1-X)*I*(1-V)*(1-Y)*Jmax))/(2*Q));
ETRabl=D*((F/2*U*(1-X)*I+V*(1-Y)*Jmax-sqrt((F/2*U*(1-X)*I+V*(1-Y)*Jmax)^2-4*Q*F/2*U*(1-X)*I*V*(1-Y)*Jmax))/(2*Q));
ETRnbl=E*((F/2*U*(1-X)*I+V*(1-Y)*Jmax-sqrt((F/2*U*(1-X)*I+V*(1-Y)*Jmax)^2-4*Q*F/2*U*(1-X)*I*V*(1-Y)*Jmax))/(2*Q));

vATPB=min(Vm_ATPB,ETRab+ETRabl)*(Bchl_ADP*Bchl_Pi-Bchl_ATP/Ke_ATPB)/(KmADP_ATPB*KmPi_ATPB*(1+Bchl_ADP/KmADP_ATPB+Bchl_Pi/KmPi_ATPB+Bchl_ATP/KmATP_ATPB+Bchl_ADP*Bchl_Pi/(KmADP_ATPB*KmPi_ATPB)));
vNADPHB=min(Vm_NADPHB,ETRnbl)*(Bchl_NADP-Bchl_NADPH/Ke_NADPHB)/(KmNADP_NADPHB*(1+Bchl_NADP/KmNADP_NADPHB+Bchl_NADPH/KmNADPH_NADPHB));
% vATPB=min(Vm_ATPB*(Bchl_ADP*Bchl_Pi-Bchl_ATP/Ke_ATPB)/(KmADP_ATPB*KmPi_ATPB*(1+Bchl_ADP/KmADP_ATPB+Bchl_Pi/KmPi_ATPB+Bchl_ATP/KmATP_ATPB+Bchl_ADP*Bchl_Pi/(KmADP_ATPB*KmPi_ATPB))),ETRab+ETRabl);
% vNADPHB=min(Vm_NADPHB*(Bchl_NADP-Bchl_NADPH/Ke_NADPHB)/(KmNADP_NADPHB*(1+Bchl_NADP/KmNADP_NADPHB+Bchl_NADPH/KmNADPH_NADPHB)),ETRnbl);
global vO2_Bchl;
vO2_Bchl=vNADPHB/2;
%Transport 17
global vOAA_M;
Km_OAA_M=0.053;
Kimal_OAA_M=7.5;
vOAA_M=Vm_OAA_M*(MC_OAA-Mchl_OAA)/(MC_OAA+Km_OAA_M*(1+MC_malate/Kimal_OAA_M));

global vMAL_M;
global vMAL;
global vMAL_B;

Km_MAL_M=0.5;
KiOAA_MAL_M=0.3;%0.3;
% Vm_MAL_B=3/20;
Km_MAL_B=1;
vMAL_M=Vm_MAL_M*(Mchl_malate-MC_malate)/(Mchl_malate+Km_MAL_M*(1+Mchl_OAA/KiOAA_MAL_M));% 2  Mchl_malate -> MC.malate
vMAL=Pmal*(MC_malate-BSC_malate);% 3 MC_malate -> BSC_malate
vMAL_B=Vm_MAL_B*(BSC_malate-Bchl_malate)/(BSC_malate+Km_MAL_B);% 4 BSC_malate -> Bchl_malate

Km_PYR_B=0.1;
vPYR_B=Vm_PYR_B*(Bchl_pyruvate-BSC_pyruvate/10)/(Km_PYR_B+Bchl_pyruvate);

global vPYR;
vPYR=Ppyr*(BSC_pyruvate-MC_pyruvate);% 6 BSC_pyruvate -> MC_pyruvate

Km_PYR_M=0.1;
vPYR_M=Vm_PYR_M*(MC_pyruvate-Mchl_pyruvate/10)/(Km_PYR_M+MC_pyruvate);

global vPEP_M;
global vPEP_B;
Km_PEP_M=0.5;
vPEP_M=Vm_PEP_M*(Mchl_PEP-MC_PEP)/(Km_PEP_M+Mchl_PEP); % 8 Mchl_PEP + MC_Pi -> MC_PEP + Mchl_Pi
Vm_PEP_B=Vm_PEP_M;
vPEP_B=Vm_PEP_B*(Bchl_PEP-BSC_PEP)/(Km_PEP_M+Bchl_PEP);

global vPGA_B;
global vDHAP_B;
global vGAP_B;
KmPGA_B = 2; 
KmGAP_B =2; 
KmDHAP_B =2; 

vPGA_B = Vtp_Bchl * Bchl_PGA/(Bchl_PGA + KmPGA_B * ( 1 + Bchl_DHAP/KmDHAP_B) * ( 1 + Bchl_GAP/KmGAP_B))-Vtp_Bchl * BSC_PGA/(BSC_PGA + KmPGA_B * ( 1 + BSC_DHAP/KmDHAP_B) * ( 1 + BSC_GAP/KmGAP_B));  
vGAP_B = Vtp_Bchl * BSC_GAP/(BSC_GAP + KmGAP_B * ( 1 + BSC_PGA/KmPGA_B) * ( 1 + BSC_DHAP/KmDHAP_B))-Vtp_Bchl * Bchl_GAP/(Bchl_GAP + KmGAP_B * ( 1 + Bchl_PGA/KmPGA_B) * ( 1 + Bchl_DHAP/KmDHAP_B));
vDHAP_B= Vtp_Bchl * BSC_DHAP/(BSC_DHAP + KmDHAP_B * ( 1 + BSC_PGA/KmPGA_B) * ( 1 + BSC_GAP/KmGAP_B))- Vtp_Bchl * Bchl_DHAP/(Bchl_DHAP + KmDHAP_B * ( 1 + Bchl_PGA/KmPGA_B) * ( 1 + Bchl_GAP/KmGAP_B));

global vPGA;
global vDHAP;
global vGAP;

vPGA=PC3P*(BSC_PGA-MC_PGA); % 10 BSC_PGA+MC.Pi -> MC.PGA+BSC_Pi
vGAP=PC3P*(MC_GAP-BSC_GAP); % 13 MC_GAP+BSC_Pi ->BSC_GAP+MC_Pi
vDHAP=PC3P*(MC_DHAP-BSC_DHAP); % 16 MC_DHAP+BSC_Pi ->BSC_DHAP+MC_Pi
global vPGA_M;
global vDHAP_M;
global vGAP_M;

KmPGA =2; 
KmGAP =2; 
KmDHAP = 2; 
vPGA_M=Vtp_Mchl * MC_PGA/(MC_PGA + KmPGA * ( 1 + MC_DHAP/KmDHAP) * ( 1 + MC_GAP/KmGAP))-Vtp_Mchl * Mchl_PGA/(Mchl_PGA + KmPGA * ( 1 + Mchl_DHAP/KmDHAP) * ( 1 + Mchl_GAP/KmGAP));
vGAP_M=Vtp_Mchl * Mchl_GAP/(Mchl_GAP + KmGAP * ( 1 + Mchl_PGA/KmPGA) * ( 1 + Mchl_DHAP/KmDHAP))-Vtp_Mchl * MC_GAP/(MC_GAP + KmGAP * ( 1 + MC_PGA/KmPGA) * ( 1 + MC_DHAP/KmDHAP));
vDHAP_M=Vtp_Mchl * Mchl_DHAP/(Mchl_DHAP + KmDHAP * ( 1 + Mchl_PGA/KmPGA) * ( 1 + Mchl_GAP/KmGAP))-Vtp_Mchl * MC_DHAP/(MC_DHAP + KmDHAP * ( 1 + MC_PGA/KmPGA) * ( 1 + MC_GAP/KmGAP));

vtATP=5*(Mchl_ATP-MC_ATP*2);

global vleakage;
global vleakage2;
vleak_B=Pco2_B*(Bchl_CO2-BSC_CO2);%0.5376*10*(Bchl_CO2-BSC_CO2);%Vleak*(Bchl_CO2-BSC_CO2);%  Bchl_CO2 -> BSC_CO2

vleak=Pco2*(BSC_CO2-MC_CO2);%Vleak*(BSC_CO2-MC_CO2);%BSC_CO2->MC_CO2
vleakage=vleak;
vleakage2=vleak_B;
% Photorespiration  
Vm_PR1=Vm_6*0.11*TempCorr_Vm_OC;%%%WY202010  
global vpr1;
vpr1=Vm_PR1*Bchl_RuBP*Bchl_O2/((Bchl_O2+KmO2_PR1*(1+Bchl_CO2/KmCO2_PR1))*(Bchl_RuBP+KmRuBP_PR1*(1+Bchl_PGA/KiPGA_PR1+Bchl_FBP/KiFBP_PR1+Bchl_SBP/KiSBP_PR1+Bchl_Pi/KiPi_PR1+Bchl_NADPH/KiNADPH_PR1)));
vpr2=Vm_PR2*Bchl_PGCA/(Bchl_PGCA+KmPGCA_PR2*(1+Bchl_GCA/KiGCA_PR2)*(1+Bchl_Pi/KiPI_PR2));
vpr3=Vm_PR3*Bper_GCA/(Bper_GCA+KmGCA_PR3);
vpr4=Vm_PR4*(Bper_GOA*Bper_GLU-Bper_KG*Bper_GLY/Ke_PS4)/((Bper_GOA+KmGOA_PS4)*(Bper_GLU+KmGLU_PS4*(1+Bper_GLY/KiGLY_PS4)));
vpr5=Vm_PR5*Bper_GLY/(Bper_GLY+KmGLY_PS5*(1+Bper_SER/KiSER_PS5));
vpr6=Vm_PR6*(Bper_GOA*Bper_SER-Bper_HPR*Bper_GLY/Ke_PR6)/((Bper_GOA+KmGOA_PR6)*(Bper_SER+KmSER_PR6*(1+Bper_GLY/KmGLY_PR6)));
vpr7=Vm_PR7*(Bper_HPR*Bper_NADH-Bper_NAD*Bper_GCEA/Ke_PR7)/(Bper_HPR+KmHPR_PR7*(1+Bper_HPR/KiHPR_PR7));
vpr9=VTgca_PR9*(Bchl_GCA/(Bchl_GCA+KmGCA_PR9*(1+Bchl_GCEA/KiGCEA_PR9))-Bper_GCA/(Bper_GCA+KmGCA_PR9*(1+Bper_GCEA/KiGCEA_PR9)));

if Pvpr8M==0
    Vm_PR8M=0;
    VTgcea_PR10M=0;
end 
if Pvpr8M==1
    Vm_PR8M=Vm_PR8;
    VTgcea_PR10M=VTgcea_PR10;
    Vm_PR8=0;
    VTgcea_PR10=0;
end 

vpr8=Vm_PR8*(Bchl_ATP*Bchl_GCEA-Bchl_ADP*Bchl_PGA/Ke_PR8)/((Bchl_ATP+KmATP_PR8*(1+Bchl_PGA/KiPGA_PR8))*(Bchl_GCEA+KmGCEA_PR8));

vpr10=VTgcea_PR10*(Bper_GCEA/(Bper_GCEA+KmGCEA_PR10*(1+Bper_GCA/KiGCA_PR10))-Bchl_GCEA/(Bchl_GCEA+KmGCEA_PR10*(1+Bchl_GCA/KiGCA_PR10)));

vpr10M=VTgcea_PR10M*(Bper_GCEA/(Bper_GCEA+KmGCEA_PR10*(1+Bper_GCA/KiGCA_PR10))-Mchl_GCEA/(Mchl_GCEA+KmGCEA_PR10));
vpr8M=Vm_PR8M*(Mchl_ATP*Mchl_GCEA-Mchl_ADP*Mchl_PGA/Ke_PR8)/((Mchl_ATP+KmATP_PR8*(1+Mchl_PGA/KiPGA_PR8))*(Mchl_GCEA+KmGCEA_PR8));


%PGA enolase, phosphoglyceromutase in MC PGA->PEP
v62=Vm_62*(MC_PGA-MC_PEP/Ke_62)/(KmPGA_62*(1+MC_PGA/KmPGA_62+MC_PEP/KmPEP_62));

KmATP_gly1=0.081;
KmF6P_gly1=0.1;
Ke_gly1=2800;
Vm_gly1B=0.035;
vgly1B=0;%Vm_gly1B*(Bchl_F6P*Bchl_ATP-Bchl_FBP*Bchl_ADP/Ke_gly1)/(KmF6P_gly1+Bchl_F6P)/(KmATP_gly1+Bchl_ATP);

%%%%%%%%%%
%PCK pathway
%%%%%%%%%%
%2.6.1.1M   PCK1
KmAsp_PCK1=2.5;KmOxog_PCK1=0.14;KmGlu_PCK1=17;KmOAA_PCK1=0.056; Ke_PCK1=1/0.148;
%2.6.1.1B   PCK2
KmAsp_PCK2=2.5;KmOxog_PCK2=0.14;KmGlu_PCK2=17;KmOAA_PCK2=0.056; Ke_PCK2=0.148;
%4.1.1.49   PCK3
KmOAA_PCK3=0.06; KmATP_PCK3=0.034*3;
%2.6.1.2B   PCK4
KmPyr_PCK4=0.33;KmGlu_PCK4=5;KmAla_PCK4=6.67;KmOxog_PCK4=0.15; Ke_PCK4=1;
%2.6.1.2M   PCK5
KmPyr_PCK5=0.33;KmGlu_PCK5=5;KmAla_PCK5=6.67;KmOxog_PCK5=0.15; Ke_PCK5=1;
%1.1.1.82B  PCK6
KmNADPH_PCK6 =0.05;  KmOAA_PCK6 =0.056;  KmNADP_PCK6 =0.045;  Kmmal_PCK6 =32.0;  Ke_PCK6 =4450.0; % No unit
global Ratio;
Vm_PCK1=1.5/20*Ratio;
Vm_PCK2=1.5/20*Ratio;
Vm_PCK4=1.5/20*Ratio;
Vm_PCK5=1.5/20*Ratio;
Vm_PCK6=0.4/20*Ratio;
Vm_PCK3=0.15/20*Ratio;
vPCK1=Vm_PCK1*(MC_OAA*MC_Glu-MC_OxoG*MC_Asp/Ke_PCK1)/(KmOAA_PCK1*KmGlu_PCK1*(1+MC_OAA/KmOAA_PCK1+MC_Glu/KmGlu_PCK1+MC_OxoG/KmOxog_PCK1+MC_Asp/KmAsp_PCK1+MC_OAA/KmOAA_PCK1*MC_Glu/KmGlu_PCK1+MC_OxoG/KmOxog_PCK1*MC_Asp/KmAsp_PCK1));
vPCK2=Vm_PCK2*(BSC_OxoG*BSC_Asp-BSC_OAA*BSC_Glu/Ke_PCK2)/(KmOxog_PCK2*KmAsp_PCK2*(1+BSC_OAA/KmOAA_PCK2+BSC_Glu/KmGlu_PCK2+BSC_OxoG/KmOxog_PCK2+BSC_Asp/KmAsp_PCK2+BSC_OAA/KmOAA_PCK2*BSC_Glu/KmGlu_PCK2+BSC_OxoG/KmOxog_PCK2*BSC_Asp/KmAsp_PCK2));
vPCK3=Vm_PCK3*BSC_OAA*BSC_ATP/((BSC_OAA+KmOAA_PCK3)*(BSC_ATP+KmATP_PCK3));
vPCK4=Vm_PCK4*(BSC_Glu*BSC_pyruvate-BSC_Ala*BSC_OxoG/Ke_PCK4)/(KmGlu_PCK4*KmPyr_PCK4*(1+BSC_Glu/KmGlu_PCK4+BSC_pyruvate/KmPyr_PCK4+BSC_Ala/KmAla_PCK4+BSC_OxoG/KmOxog_PCK4+BSC_Glu/KmGlu_PCK4*BSC_pyruvate/KmPyr_PCK4+BSC_Ala/KmAla_PCK4*BSC_OxoG/KmOxog_PCK4));
vPCK5=Vm_PCK5*(MC_Ala*MC_OxoG-MC_Glu*MC_pyruvate/Ke_PCK5)/(KmAla_PCK5*KmOxog_PCK5*(1+MC_Glu/KmGlu_PCK5+MC_pyruvate/KmPyr_PCK5+MC_Ala/KmAla_PCK5+MC_OxoG/KmOxog_PCK5+MC_Glu/KmGlu_PCK5*MC_pyruvate/KmPyr_PCK5+MC_Ala/KmAla_PCK5*MC_OxoG/KmOxog_PCK5));
vPCK6=Vm_PCK6*(Bchl_OAA*Bchl_NADPH-Bchl_NADP*Bchl_malate/Ke_PCK6)/( KmOAA_PCK6* KmNADPH_PCK6*(1+Bchl_OAA/KmOAA_PCK6+ Bchl_NADPH/KmNADPH_PCK6+ Bchl_NADP/KmNADP_PCK6+ Bchl_malate/Kmmal_PCK6+ Bchl_OAA*Bchl_NADPH/(KmOAA_PCK6* KmNADPH_PCK6)+ Bchl_NADP*Bchl_malate/(KmNADP_PCK6* Kmmal_PCK6)));
%%%%%%%%%%%%%
%Transport for PCK
%%%%%%%%%%%%%
%Tasp
global vAsp;
vAsp=0.664/20*(MC_Asp-BSC_Asp);
%TAla
vAla=0.8715/20*(BSC_Ala-MC_Ala);
%TPEP
vPEP=0.653625/20*(BSC_PEP-MC_PEP);
%TOAAB
vOAA_B=1*(BSC_OAA-Bchl_OAA);
%ATPB
VtATP=3;
vATP_B=VtATP*(Bchl_ATP-2*BSC_ATP);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NAD-ME pathway
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.1.1.37 MDH_BM
%KeMDH_BM=1/(6.65*10^(-5));
Vm_MDH_BM=0.3;
KmOAA_MDH_BM=0.126;%Maize
KmNADH_MDH_BM=0.025;
%1.1.1.38 NADME
VmNADME=2;
%KeNADME=0.0011*55.35*1000;% Wedding, R.T.; Black, M.K.; Plant Physiol.; 72, 1021 (1983).
Km_MAL_NADME=3;
Km_NAD_NADME=0.5;
Ki_NADH_NADME=0.15;
Ki_PYR_NADME=14;

vMDH_BM=Vm_MDH_BM*Bmito_OAA*Bmito_NADH/(KmOAA_MDH_BM+Bmito_OAA)/(KmNADH_MDH_BM+Bmito_NADH);
vNADME=VmNADME*Bmito_MAL*Bmito_NAD/(Km_MAL_NADME*(1+Bmito_PYR/Ki_PYR_NADME)+Bmito_MAL)/(Km_NAD_NADME*(1+Bmito_NADH/Ki_NADH_NADME)+Bmito_NAD);
vtOAA_Bm=1*(BSC_OAA-Bmito_OAA);
vPYR_Bm=1*(Bmito_PYR-BSC_pyruvate);
vMAL_Bm=1*(BSC_malate-Bmito_MAL);

%%%%%O2 diffusion%%%%
global vtO2;
global vtO2_B;
global vtO2_M;

vtO2=0.122633*(BSC_O2-MC_O2);
vtO2_B=Pco2_B*(Bchl_O2-BSC_O2);
vtO2_M=Pco2_B*(Mchl_O2-MC_O2);
Bchl_Asp=BSC_Asp;
Bchl_Ala=BSC_Ala;
Mchl_Asp=MC_Asp;
Mchl_Ala=MC_Ala;

vAspB=10*(BSC_Asp-Bchl_Asp);
vAlaB=10*(BSC_Ala-Bchl_Ala);
vAspM=10*(MC_Asp-Mchl_Asp);
vAlaM=10*(MC_Ala-Mchl_Ala);

global pathway_option;
if pathway_option==0
    vPCK1=0;
    vPCK2=0;
    vPCK3=0;
    vPCK4=0;
    vPCK5=0;
    vPCK6=0;
    vAsp=0;
    vAla=0;
%    vPEP=0;
    vOAA_B=0;
    vATP_B=0;
    vMDH_BM=0;
    vNADME=0;
    vtOAA_Bm=0;
    vPYR_Bm=0;
    vMAL_Bm=0;
end
if pathway_option==1
   vPCK3=0;
   vATP_B=0;
   vPEP=0;
   vMDH_BM=0;
   vNADME=0;
   vtOAA_Bm=0;
   vPYR_Bm=0;
   vMAL_Bm=0;
end
if pathway_option==2
   vPCK6=0; 
   vOAA_B=0;
   vMDH_BM=0;
   vNADME=0;
   vtOAA_Bm=0;
   vPYR_Bm=0;
   vMAL_Bm=0;
end
if pathway_option==3
    vMDH_BM=0;
    vNADME=0;
    vtOAA_Bm=0;
    vPYR_Bm=0;
    vMAL_Bm=0;
end
if pathway_option==4
   v3=0;
   v4=0;
   v5=0;
   vPCK6=0;
   vPYR_B=0;
   vPYR_M=0;
   vOAA_B=0;
   vOAA_M=0;
   vMAL_M=0;
   vMAL_B=0;
   vPEP_M=0;
   vMDH_BM=0;
   vNADME=0;
   vtOAA_Bm=0;
   vPYR_Bm=0;
   vMAL_Bm=0;
end

if pathway_option==6
    vMAL_B=0;
    vMDH_BM=0;
    vNADME=0;
    vtOAA_Bm=0;
    vPYR_Bm=0;
    vMAL_Bm=0;
end

if pathway_option==7
    v3=0;
    vOAA_M=0;
    vMAL_M=0;
    vMAL=0;
    vMAL_B=0;
    vPYR_B=0;
    v4=0;
    vPCK6=0;
    vOAA_B=0;
    vPCK3=0;
    vMAL_Bm=0;
    vPYR=0;
end

if pathway_option==8
    vtOAA_Bm=0;
    vMDH_BM=0;
    vOAA_B=0;
    v4=0;
    vPCK6=0;
    vMAL_B=0;
    vPYR_B=0;
    vNAE=vNADME;% 1NAD-ME-> 2.5 ATP
end

if Para_mata==0 % C4 collatz 1992
Q10Temperature =PhotosynthesQ10^((LeafTemperature - 25.0) / 10.0);
Rd = Rd25 * Q10Temperature / (1.0 + exp(1.3 * (LeafTemperature - 55.0)));
ThetaC4=0.83;%curvature parameter
BetaC4=0.93;%curvature parameter
kC4=0.7;% initial slope of photosynthetic CO2 response
kiC4=kC4*2^((LeafTemperature-25)/10);
%VmaxC4=45;%umol m-2 s-1 maximum rubisco capacity
VmaxiC4=VmaxC4*2^((LeafTemperature-25)/10)/((1+exp(0.3*(13-LeafTemperature)))*(1+exp(0.3*(LeafTemperature-36))));
Ji=0.05*Convert * Radiation_PAR;%alpha*ar*f*Qp;
Jc=Ci*10^-6*Pressure*kiC4*10^6/Pressure;%pi*(kp-L/pi)/P;
Je=VmaxiC4;
%Ax=[Ji Jc Je];
M=(Je+Ji-sqrt((Je+Ji)^2-4*ThetaC4*Je*Ji))/(2*ThetaC4);
GrossAssimilation=(M+Jc-sqrt((M+Jc)^2-4*BetaC4*M*Jc))/(2*BetaC4);
NetAssimilation=GrossAssimilation-Rd;
%end
end
if Para_mata==1 % metabolic model
NetAssimilation=vinf*1000;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 LeafDimension = 0.028; %Leaf width/needle diameter m
 CForced = 4.322 / 1000.0; 
 CFree = 1.6361 / 1000.0;
 TemperatureKelvin = WeatherTemperature + 273.15;% Air temperature K
 LeafTemperatureKelvin = LeafTemperature + 273.15;%Leaf temperature K
 ESatweather = 0.611 * exp(17.502 * WeatherTemperature / (WeatherTemperature + 240.97));%Vapor pressure kPa
 Ea = WeatherRH * ESatweather; %Vapor pressure kPa
 ESaturation = 0.611 * exp(17.502 * LeafTemperature / (LeafTemperature + 240.97));%Vapor pressure kPa

 %Forced convection
 GbForced = CForced * TemperatureKelvin^0.56 * sqrt((TemperatureKelvin + 120.0)* (WeatherWind / LeafDimension / Pressure)); %// m/s
 %Free convection
 %GbFree = GbForced;% m/s
 %Eb = (LeafGs / 41.1 * Ei + GbFree * Ea*1000) / (LeafGs / 41.1 + GbFree); % Stomatal conductance from moles/m2 leaf area/s to m/s
 TDifference = (LeafTemperatureKelvin / (1.0 - 0.378 * Eb /Pressure)) -(TemperatureKelvin / (1.0 - 0.378 * Ea*1000 / Pressure));
 GbFree = CFree * LeafTemperatureKelvin^0.56 * sqrt((LeafTemperatureKelvin + 120.0) /Pressure) * (abs(TDifference) / LeafDimension)^0.25;% m/s

% Maximum of two conductances
 if GbFree >= GbForced
     Gbw = GbFree;
 else
     Gbw = GbForced; % m/s
 end
 Gbw = Gbw * 41.4; %Conversion from m/s to moles/m2 leaf area/s
 Gb=Gbw/1.37;%or 1.37
 Cb = Air_CO2 - 1.37 * NetAssimilation / Gb; % ppm ***********
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Update stomatal conductance moles/m2 leaf area/s
Gsw0 = BallBerryInterceptC4 + BallBerrySlopeC4*NetAssimilation * Eb / ESaturation/Cb;  
a = -0.1081;%Jarvist 1976
b = 1.009;
c = 1.104;
WaterStressFactor=a*exp(-b*PhiLeaf)+c;

Gsw0 = gsxsen*Gsw0 * WaterStressFactor; %Apply water stress factor
Gs0=Gsw0/1.6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if gs response time is 0
if Gs0<0
    Gs0=0;
end

%if gs response time is not O
if GsResponse==1
    if kdcon==1
      kd_Gs= 1/(kd/60);
      ki_Gs= 1/(ki/60);
    end
    if kdcon==0
      kd_Gs= 1/((kd+ainter*I)/60);
      ki_Gs= 1/(ki/60);
    end
      
if Gs<Gs0
vgs=(Gs0-Gs)/ki_Gs;
end
if Gs>=Gs0
vgs=(Gs0-Gs)/kd_Gs;
end
end
if GsResponse==0
Gs=Gs0;
Gsw=Gsw0;
vgs=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gctotal = Gb * Gs / (Gb+Gs); % moles/m2 leaf area/s
Gwtotal = Gbw * Gsw / (Gbw+Gsw);
Cb= Air_CO2 - NetAssimilation / Gb;
%Eb= NetAssimilation*(Pressure / 1000.0)/ (10^6.0*Gbw)+Ea;
Eb=Gwtotal*(ESaturation-Ea)/Gbw+Ea;
vCO2b=Gb*(Air_CO2-Cb);
vCO2s=Gs*(Cb-Ci);
vCO2total=Gctotal*(Air_CO2-Ci);
vH2Ob=Gbw*(Eb-Ea)/(Pressure / 1000.0)*10^6.0;
vH2Os=Gsw*(ESaturation-Eb)/(Pressure / 1000.0)*10^6.0;
vH2Ototal=Gwtotal*(ESaturation-Ea)/(Pressure / 1000.0)*10^6.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%leaf energy balance
 Epsilon = 0.96;% Leaf thermal emissivity
 LEFactor = 1.0; %%%%1.0%%%%%
 LWFactor = 2.0; 
 HFactor = 2.0;%%%%%GGGGGGGG
 LeafEnergyFluxMe = 0.506 * NetAssimilation; %Energy in biochemical reactions W/m2
 if Radiation_LW == 0
 Radiation_LW = LWFactor * Epsilon * Boltzman *(273.15 + WeatherTemperature) ^4;
 end
 SensibleHeat = HFactor * ConstantsCp * 0.924 * Gbw * (LeafTemperature - WeatherTemperature);

 Emission = LWFactor * Epsilon * Boltzman * (273.15 + LeafTemperature)^4.0;
 LatentHeat = LEFactor *vH2Ototal/10^6.0* LatentHeatVaporization;
 Abs=0.85;
 EnergyBalanceResidual = Radiation_PAR*Abs + Radiation_NIR+Radiation_LW -Emission -SensibleHeat -LatentHeat - LeafEnergyFluxMe;


global reaction_flux;
Reaction_v=zeros(9,1);
Reaction_v(1)=NetAssimilation;
Reaction_v(2)=vCO2b;
Reaction_v(3)=vCO2s;
Reaction_v(4)=vH2Ob;
Reaction_v(5)=vH2Os;
Reaction_v(6)=EnergyBalanceResidual;
Reaction_v(7)=vCO2total;
Reaction_v(8)=vH2Ototal;
Reaction_v(9)=vgs;
reaction_flux=Reaction_v;


if (TIME_N ==0)
    TIME_N = 1;
end

% % %if (t > OLD_TIME)
% %     TIME_N = TIME_N + 1;
% %     OLD_TIME = t;%Gs_VEL(TIME_N-2,1) ;
% % %end


if (t > OLD_TIME)
    TIME_N = TIME_N + 1;
    OLD_TIME = t;%Gs_VEL(TIME_N-2,1) ;
end
Gs_VEL(TIME_N,1) = t;
Gs_VEL(TIME_N,2) = NetAssimilation;
Gs_VEL(TIME_N,3) = vCO2b;
Gs_VEL(TIME_N,4) = vCO2s;
Gs_VEL(TIME_N,5) = vH2Ob;
Gs_VEL(TIME_N,6) = vH2Os;
Gs_VEL(TIME_N,7) = EnergyBalanceResidual;
Gs_VEL(TIME_N,8) = vCO2total;
Gs_VEL(TIME_N,9) = vH2Ototal;
Gs_VEL(TIME_N,10) = vgs;
Gs_VEL(TIME_N,11) = vleak/v2;
Gs_VEL(TIME_N,12) =Gs;
Gs_VEL(TIME_N,13) =Cb;
Gs_VEL(TIME_N,14) =Ci;
Gs_VEL(TIME_N,15) =Gbw;
global enzyme_flux;
Enz_v=zeros(102,1);   
Enz_v(1)=v1;
Enz_v(2)=v2;
Enz_v(3)=v3;
Enz_v(4)=v4;
Enz_v(5)=v5;
Enz_v(6)=v6;
Enz_v(7)=v7;
Enz_v(8)=v8;
Enz_v(9)=v10;
Enz_v(10)=v11;
Enz_v(11)=v12;
Enz_v(12)=v13;
Enz_v(13)=v14;
Enz_v(14)=v15;
Enz_v(15)=v18;
Enz_v(16)=v7Mchl;
Enz_v(17)=v8Mchl;
Enz_v(18)=vStarch1;
Enz_v(19)=vPGASink;
Enz_v(20)=vSuc1;
Enz_v(21)=vSuc2;
Enz_v(22)=vSuc3;
Enz_v(23)=vSuc4;
Enz_v(24)=vSuc7;
Enz_v(25)=vSuc8;
Enz_v(26)=vSuc9;
Enz_v(27)=vSuc10;
Enz_v(28)=vATPM;
Enz_v(29)=vNADPHM;
Enz_v(30)=vATPB;
Enz_v(31)=vOAA_M;
Enz_v(32)=vMAL_M;
Enz_v(33)=vMAL;
Enz_v(34)=vMAL_B;
Enz_v(35)=vPYR_B;
Enz_v(36)=vPYR;
Enz_v(37)=vPYR_M;
Enz_v(38)=vPEP_M;
Enz_v(39)=vPGA_B;
Enz_v(40)=vPGA;
Enz_v(41)=vPGA_M;
Enz_v(42)=vGAP_M;
Enz_v(43)=vGAP;
Enz_v(44)=vGAP_B;
Enz_v(45)=vDHAP_M;
Enz_v(46)=vDHAP;
Enz_v(47)=vDHAP_B;
Enz_v(48)=vleak_B;
Enz_v(49)=vleak;
Enz_v(50)=vNADPHB;
Enz_v(51)=vpr1;
Enz_v(52)=vpr2;
Enz_v(53)=vpr3;
Enz_v(54)=vpr4;
Enz_v(55)=vpr5;
Enz_v(56)=vpr6;
Enz_v(57)=vpr7;
Enz_v(58)=vpr8;
Enz_v(59)=vpr9;
Enz_v(60)=vpr10;
Enz_v(61)=vinf;
Enz_v(62)=v62;
Enz_v(63)=v78;
Enz_v(64)=v78Mchl;
Enz_v(65)=vStarch2;
Enz_v(66)=vgly1B;
Enz_v(67)=vhexp;
Enz_v(68)=vSta1;
Enz_v(69)=vSta2;
Enz_v(70)=vSta3;
Enz_v(71)=vtATP;
Enz_v(72)=vPCK1;
Enz_v(73)=vPCK2;
Enz_v(74)=vPCK3;
Enz_v(75)=vPCK4;
Enz_v(76)=vPCK5;
Enz_v(77)=vPCK6;
Enz_v(78)=vAsp;
Enz_v(79)=vAla;
Enz_v(80)=vPEP;
Enz_v(81)=vOAA_B;
Enz_v(82)=vATP_B;
Enz_v(83)=vO2_Mchl;
Enz_v(84)=vO2_Bchl;
Enz_v(85)=vtO2;
Enz_v(86)=vtO2_B;
Enz_v(87)=vtO2_M;
Enz_v(88)=v5B;
Enz_v(89)=vPEP_B;
Enz_v(90)=vpr10M;
Enz_v(91)=vpr8M;
Enz_v(92)=vMDH_BM;
Enz_v(93)=vNADME;
Enz_v(94)=vtOAA_Bm;
Enz_v(95)=vPYR_Bm;
Enz_v(96)=vMAL_Bm;
Enz_v(97)=vAspB;
Enz_v(98)=vAlaB;
Enz_v(99)=vAspM;
Enz_v(100)=vAlaM;
Enz_v(101)=vEA_PPDKRP_I;
Enz_v(102)=vEA_PPDKRP_A;
enzyme_flux=Enz_v;

EnzAct_v=zeros(11,1);
EnzAct_v(1)=vATPsynthase_Act_Mchl;
EnzAct_v(2)=vNADPMDH_Act;
EnzAct_v(3)=vGAPDH_Act_Mchl;
EnzAct_v(4)=vATPsynthase_Act_Bchl;
EnzAct_v(5)=vPEPC_Act;
EnzAct_v(6)=ActRubisco0;
EnzAct_v(7)=vRubisco_Act;
EnzAct_v(8)=vGAPDH_Act_Bchl;
EnzAct_v(9)=vFBPase_Act;
EnzAct_v(10)=vSBPase_Act;
EnzAct_v(11)=vPRK_Act;
EnzAct_v(12)=vRCA_Act;

LMEnz_v(1:9,1)=Reaction_v;
LMEnz_v(10:111,1)=Enz_v;
LMEnz_v(112:123,1)=EnzAct_v;


if (TIME_M ==0)
    TIME_M = 1;
end

if (t > OLD_TIME_M)
    TIME_M = TIME_M + 1;
    OLD_TIME_M = t;
end
Meta_VEL(TIME_M,1) = t;
Meta_VEL(TIME_M,2) = vinf;
Meta_VEL(TIME_M,3) =v1;
Meta_VEL(TIME_M,4) =v2;
Meta_VEL(TIME_M,5) =v3;
Meta_VEL(TIME_M,6) =v4;
Meta_VEL(TIME_M,7) =v5;
Meta_VEL(TIME_M,8) =v6;
Meta_VEL(TIME_M,9) =vpr1;
end


