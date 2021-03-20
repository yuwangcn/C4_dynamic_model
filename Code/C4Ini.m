%the code for the initalizer EnzymeIni.m
%simulation of our enzyme system
 
function con = C4Ini()%(Begin, GenNum)
% global I;
global VPPDKsen;
global VMDHsen;
global VMEsen;
global VGAPDHsen;
global VSBPsen;
global VFBPsen;
global VPRKsen;
global Vexsen;
global Jmaxsen;

Convert=1E6/(2.35E5); 
global Bchl_CA;
global Bchl_CN;
Bchl_CA= 1.5;
Bchl_CN= 0.5;
%Bchl_CP= 40.0;
global MC_CU;
global MC_CA;
MC_CU=1.5;
MC_CA=1.0;
global Mchl_CA;
global Mchl_CN;
Mchl_CA= 1.5;
Mchl_CN= 0.5;

global MC_UTP;
MC_UTP=1.26;%0.75

global Bper_GLU;
global Bper_KG;
global Bper_NADH;
global Bper_NAD;
Bper_GLU=24;
Bper_KG=0.4;
Bper_NADH=0.47;
Bper_NAD=0.4;

%MC

MC_HCO3= 0.3;%WY1911 0.005
MC_OAA=0.01;
MC_PEP=0.1;
MC_Malate=1.0;
MC_Pyruvate=2;%%%%%%0.02,0427%%%%%%
MC_PGA=0.3;
MC_FBP=0.04;
MC_UDPG=0.035;
MC_SUCP=0.0;
MC_SUC=0.0;
MC_F26BP=7.8E-5;
MC_ATP=0.35;%0.35
MC_T3P=0.55;
MC_HexP=2.4;
MC_Sucrose=0.0;
%Mchl_
Mchl_OAA= 0.005;
Mchl_Malate =1.8;
Mchl_PEP =0.1;
Mchl_Pyruvate= 0.01;
Mchl_NADPH= 0.21;
Mchl_ATP= 1.4;
Mchl_PGA= 0.04;
Mchl_DPGA= 0.0;
Mchl_T3P= 0.6;
%BSC
BSC_T3P= 0.45;
BSC_PGA= 0.2;
BSC_Malate=0.8;
BSC_Pyruvate= 2;%%%0427,0.15%%%%
BSC_CO2= 0.001;
%Bchl
Bchl_CO2= 0.01;
Bchl_RuBP= 2.0;
Bchl_PGA= 0.3;
Bchl_DPGA= 0;
Bchl_ATP= 1.4;
Bchl_NADPH= 0.1;
Bchl_SBP= 0.015;
Bchl_S7P= 0.045;
Bchl_FBP= 0.06;
Bchl_E4P= 0.05;
Bchl_Starch= 0.0;
Bchl_Rubisco= 1.456965457;
Bchl_T3P= 0.5;
Bchl_HexP= 6;%2.2;WY2002
Bchl_Pent =0.05;
Bchl_Malate= 0.3;
Bchl_Pyruvate= 0.23;

Bchl_PGCA=0.0029;
Bchl_GCA=0.36;
Bchl_GCEA=0.1812;

Bper_GCA=0.36;
Bper_GOA=0.028;
Bper_GLY=1.8;
Bper_SER=7.5;
Bper_HPR=0.0035;
Bper_GCEA=0.1812;

Bchl_PPi=0;
Bchl_ADPG=0;


global CI;
MC_CO2=0.8*CI;

MC_Glu=15;
MC_OxoG=3;
MC_Asp=5*2;
MC_Ala=5*2;

BSC_OxoG=3;
BSC_Glu=15;
BSC_Asp=5*2;
BSC_Ala=5*2;
BSC_OAA=0;
BSC_PEP=0.1;
BSC_ATP=0.5;
Bchl_OAA=0;

global O2;
MC_O2=O2;
Mchl_O2=O2;
BSC_O2=O2;
Bchl_O2=O2;

Bchl_PEP=0.1;
Mchl_GCEA=0.1812;

Bmito_OAA=0;
Bmito_MAL=0;
Bmito_PYR=5;
Bmito_CO2=0.001;
Bmito_NADH=0.3;

Bchl_Asp=0;
Bchl_Ala=0;
Mchl_Asp=0;
Mchl_Ala=0;
E_PPDK_Mchl=0;% E_PPDK active
EP_PPDK_Mchl=0.616;%1.28;%EP_PPDK inactive %WY202010 0.616 asumming Vmax_PPDK = 80


con=zeros(1,81);
con(1)=MC_HCO3;
con(2)=MC_OAA;
con(3)=MC_PEP;
con(4)=MC_Malate;
con(5)=MC_Pyruvate;
con(6)=MC_PGA;
con(7)=MC_FBP;
con(8)=MC_UDPG;
con(9)=MC_SUCP;
con(10)=MC_SUC;
con(11)=MC_F26BP;
con(12)=MC_ATP;
con(13)=MC_T3P;
con(14)=MC_HexP;
con(15)=MC_Sucrose;
con(16)=Mchl_OAA;
con(17)=Mchl_Malate;
con(18)=Mchl_PEP;
con(19)=Mchl_Pyruvate;
con(20)=Mchl_NADPH;
con(21)=Mchl_ATP;
con(22)=Mchl_PGA;
con(23)=Mchl_DPGA;
con(24)=Mchl_T3P;
con(25)=BSC_T3P;
con(26)=BSC_PGA;
con(27)=BSC_Malate;
con(28)=BSC_Pyruvate;
con(29)=BSC_CO2;
con(30)=Bchl_CO2;
con(31)=Bchl_RuBP;
con(32)=Bchl_PGA;
con(33)=Bchl_DPGA;
con(34)=Bchl_ATP;
con(35)=Bchl_NADPH;
con(36)=Bchl_SBP;
con(37)=Bchl_S7P;
con(38)=Bchl_FBP;
con(39)=Bchl_E4P;
con(40)=Bchl_Starch;
con(41)=Bchl_Rubisco;
con(42)=Bchl_T3P;
con(43)=Bchl_HexP;
con(44)=Bchl_Pent;
con(45)=Bchl_Malate;
con(46)=Bchl_Pyruvate;

con(47)=Bchl_PGCA;
con(48)=Bchl_GCA;
con(49)=Bchl_GCEA;

con(50)=Bper_GCA;
con(51)=Bper_GOA;
con(52)=Bper_GLY;
con(53)=Bper_SER;
con(54)=Bper_HPR;
con(55)=Bper_GCEA;
con(56)=MC_CO2;

con(57)=Bchl_PPi;
con(58)=Bchl_ADPG;

con(59)=MC_Glu;
con(60)=MC_OxoG;
con(61)=MC_Asp;
con(62)=MC_Ala;
con(63)=BSC_OxoG;
con(64)=BSC_Glu;
con(65)=BSC_Asp;
con(66)=BSC_Ala;
con(67)=BSC_OAA;
con(68)=BSC_PEP;
con(69)=BSC_ATP;
con(70)=Bchl_OAA;
con(71)=MC_O2;
con(72)=Mchl_O2;
con(73)=BSC_O2;
con(74)=Bchl_O2;
con(75)=Bchl_PEP;%%%%%%%WY PPDK in BSCytosol
con(76)=Mchl_GCEA;

con(77)=Bmito_OAA;
con(78)=Bmito_MAL;
con(79)=Bmito_PYR;
con(80)=Bmito_CO2;
con(81)=Bmito_NADH;

con(82)=Bchl_Asp;
con(83)=Bchl_Ala;
con(84)=Mchl_Asp;
con(85)=Mchl_Ala;
con(86)=E_PPDK_Mchl;% E_PPDK active
con(87)=EP_PPDK_Mchl;%EP_PPDK inactive

 
%Km Ki(mM)
%4.2.1.1        1
KmCO2_1 = 2.8;%2.8;  
Ke_1 =11.2;%20;%;% No unit  k=4.3*10^(-7)  PH=7.5 13.6 11.82
%4.1.1.31       2 %KmHCO3_2 =0.05 KmPEP_2 =0.047
% 
KmPEP_2 =0.1; KmHCO3_2 =0.02; Kimal_2 =1;
%%%%%WY202003 Maize KmHCO3_2 =0.04216; KmPEP_2 =1;   Kimal_2 =15; 
%KmPEP_2 =0.038; 
%1.1.1.82       3
KmNADPH_3 =0.024;  KmOAA_3 =0.056;  KmNADP_3 =0.073;  Kmmal_3 =32.0;  
Ke_3 =4450.0; % No unit
%1.1.1.40       4 Kmmal_4 =0.04;
KmCO2_4 =1.1;  KmNADP_4 =0.0080;  KmNADPH_4 =0.045;  KmPyr_4 =3;  Kmmal_4 =0.23;  Ke_4 =0.051*55.35*1000;%0.0344 KmNADP_4 =0.008
%2.7.9.1        5 KiPEP_5 =0.5;
KiPEP_5 =0.15;  KmATP_5 =0.082;  KmPyr_5 =0.082;% KmPyr_5 =0.25
%4.1.1.39       6  KiPGA_6 =0.84; KiFBP_6 =0.4;
KmCO2_6 =0.0162;  KmO2_6 =0.183;  KmRuBP_6 =0.02;  KiPGA_6 =2.52;  KiFBP_6 =0.04;  KiSBP_6 =0.075;  KiPi_6 =0.9*3;  KiNADPH_6 =0.07*3;%KiNADPH_6 =0.07;KiPi_6 =0.9
%2.7.2.3        7
KmADP_7 =0.5;  KmATP_7 =0.3;  KmPGA_7 =2.4;%1.2;%KmPGA_7 =0.24
%1.2.1.13       8
KmDPGA_8 =0.4;  KmNADPH_8=0.1;
%5.3.1.1        9
Ke_9=0.05;    
%4.1.2.13FBP    10
KmDHAP_10 =0.45;  KmFBP_10 =0.0923;  KmGAP_10 =0.04;  %KmGAP_10 =0.3; 
Ke_10 =7.1; % 1/millimolarity
%3.1.3.11       11 KmFBP_11 =0.033; Chloroplast fructose-1,6-bisphosphatase
%from Oryza differs in salt tolerance property from the Porteresia enzyme
%and is protected by osmolytes  Ghosh, S.; Bagchi, S.; Lahiri Majumder, A.;Plant Sci. 160, 1171-1181 (2001)
KiF6P_11 =0.7;  KiPi_11 =12.0;  KmFBP_11 =0.066;  Ke_11 =666000.0;
%4.1.2.13SBP    12
KmDHAP_12 =0.4;  KmE4P_12 =0.2;  
Ke_12 =1.017;  % 1/millimolarity
%3.1.3.37       13
KiPi_13 =12.0;  KmSBP_13 =0.05;  Ke_13 =666000.0;
%2.2.1.1X       14 Datta, A.G.; Racker, E.; J. Biol. Chem.; 236, 617 (1961). 
KmE4P_14 =0.1;  KmF6P_14 =0.1;  KmGAP_14 =0.1;  KmXu5P =0.1;  
Ke_14 =0.084;  % No unit Ke_14 =10.0£»Ke_14 =0.076
%2.2.1.1R       15
KmGAP_15 =0.072;  KmRi5P_15 =1.5;  KmS7P_15 =0.015;  KmXu5P_15 =0.1;  % KmS7P_15 =0.015  KmS7P_15 =0.46
Ke_15 =0.9;%1.176470588235294;  % No unit Datta, A.G.; Racker, E.; J. Biol. Chem.; 236, 617 (1961). 
%5.3.1.6:Chl    16
Ke_16=0.4;
%5.1.3.1:Chl    17
Ke_17=0.67;
%2.7.1.19       18 % KiPGA_18 =1.0;  KiPi_18 =1.0; KmATP_18 =0.059;KmATP_18=0.625
KiADP_18 =2.5;  Ki_ADP_18 =0.4;  KiPGA_18 =2.0;  KiPi_18 =4.0;  KiRuBP_18 =0.7;  KmATP_18 =0.625;  KmRu5P_18 =0.05; % KiPi_18 =4.0; 
Ke_18 =6846.0; % No unit

%2.7.2.3:MChl  7Mchl
KmADP_7Mchl =0.5;  KmATP_7Mchl =0.3;  KmPGA_7Mchl =2.4;
%1.2.1.13:MChl   8Mchl
KmDPGA_8Mchl =0.4;  KmNADPH_8Mchl =0.1;  

%StarchSynthesis:Chl   Starchsyn  KmG1P_Starch =0.24;  KaF6P_Starch =0.06;
%KaFBP_Starch =0.06;
KiADP_Starch =10.0;  KmATP_Starch =0.08;  KmG1P_Starch =0.48;  KaF6P_Starch =0.12;  KaFBP_Starch =0.12; KaPGA_Starch =0.3;% Ka (No unit)KmG1P_Starch =0.08;  KaF6P_Starch =0.02;  KaFBP_Starch =0.02; KaPGA_Starch =0.1
Ke_Starch1=2.3;   Ke_Starch2=0.058;
%PGASink
KmPGA_PGASink=1;%0.5;

%4.1.2.13FBP:Cel    Suc1
KmDHAP_Suc1 =0.45;  KmGAP_Suc1 =0.04/2;  KmFBP_Suc1 =0.0023;   %KmDHAP_Suc1 =0.4;  KmGAP_Suc1 =0.1;  KmFBP_Suc1 =0.2;
Ke_Suc1 =12;%12.0; 10000 % 1/millimolarity
%3.1.3.11:Cel    Suc2 KiF26BP_Suc2 =7.0E-5;KmFBP_Suc2 =0.165;Ke_Suc2=6663.0
%Roles of the residues Lys115 and Tyr116 in the binding of an allosteric inhibitor AMP to pea cytosolic D-fructose-1,6-bisphosphatase
%Jang, H.; Cho, M.; Kwon, Y.; Bhoo, S.H.; Jeon, J.; Hahn, T.; J. Appl. Biol. Chem. 51, 45-49 (2008)
KiF26BP_Suc2 =0.00007;  KiF6P_Suc2 =0.7;  KiPi_Suc2 =12.0;  KmFBP_Suc2 =0.00108;  Ke_Suc2 =174.0;
%5.3.1.9:Cel    Suc5
Ke_Suc5=2.3;
%5.4.2.2:Cel    Suc6
Ke_Suc6=0.0584;
%2.7.7.9:Cel    Suc7
KmG1P_Suc7 =0.14;  KmPPi_Suc7 =0.11;  KmUDPG_Suc7 =0.12;  KmUTP_Suc7 =0.1;  
Ke_Suc7 =0.31;%0.31; % No unit
%2.4.1.14:Cel   Suc8
KiFBP_Suc8 =0.8;  KiPi_Suc8 =5.0;  KiSuc_Suc8 =50.0;  KiSucP_Suc8 =0.4;  KiUDP_Suc8 =0.7;  KmF6P_Suc8 =0.8;  KmUDPG_Suc8 =1.3;  
Ke_Suc8 =10.0;  % No unit
%3.1.3.24:Cel     Suc9
KmSuc_Suc9 =80.0;  KmSucP_Suc9 =0.35;  Ke_Suc9 =780.0;
%SUCSink:Cel   Suc10
KmSuc_Suc10 =1.5;
%2.7.1.105:Cel     Suc3
KiADP_Suc3 =0.16;  KIDHAP_Suc3 =0.7;  KmATP_Suc3 =1.32;  KmF26BP_Suc3 =0.021;  KmF6P_Suc3 =1.4;  
Ke_Suc3 =590.0;% No unit
%3.1.3.46:Cel    Suc4
KiF6P_Suc4 =0.1;  KiPi_Suc4 =0.5*10;  KmF26BP_Suc4= 0.032; 
%3.6.1.1:Cel 
KePi=128.4;
%3.6.3.14:MChl   ATPM
KmADP_ATPM =0.014;  KmATP_ATPM =0.11;  KmPi_ATPM =0.3;
Ke_ATPM =5.734;        %1/millimolarity
global Xd;%0.667
X =0.667;  Y =0.6;  F =0.7225;  Q =0.9;  D =1;% No unit
%1.18.1.2:MChl  NADPHM
KmNADP_NADPHM =0.05;  KmNADPH_NADPHM =0.058; % KmNADP_NADPHM =0.0072;  KmNADPH_NADPHM =0.0058;  
Ke_NADPHM =502;  E =0.5; % No unit
%V3.6.3.14:Chl       ATPB
KmADP_ATPB =0.014;  KmPi_ATPB =0.11;   KmATP_ATPB =0.3;
Ke_ATPB =5.734; % 1/millimolarity
G =0.667; % No unit
%1.18.1.2:BChl  NADPHB
KmNADP_NADPHB =0.05;  KmNADPH_NADPHB =0.058;  
Ke_NADPHB =502;  % No unit

%4.1.1.39 O2 1
%KmCO2_6 =0.01935;  KmO2_6 =0.222;  KmRuBP_6 =0.02;  KiPGA_6 =2.52; KiFBP_6 =0.8;  KiSBP_6 =0.75;  KiPi_6 =0.9;  KiNADPH_6 =0.07;
KmCO2_PR1=0.0162;  KmO2_PR1=0.183;  KmRuBP_PR1=0.02;  KiPGA_PR1=2.52;  KiFBP_PR1=0.04;  KiSBP_PR1=0.75;  KiPi_PR1=0.9*3;  KiNADPH_PR1=0.21;
%3.1.3.18 2
KmPGCA_PR2=0.026;  KiPI_PR2=2.55;  KiGCA_PR2=94.0;%KmPGCA_PR2=0.026; 0.57
%1.1.3.15 3
KmGCA_PR3= 0.1;%KmGCA_PR3= 0.1;0.02
%2.6.1.4 4
Ke_PS4= 607.0;  KmGOA_PS4=0.15;  KmGLU_PS4= 1.7;  KiGLY_PS4=2.0;
%GLY_SER:Mit 5
KmGLY_PS5= 6.0;  KiSER_PS5=4.0;
%2.6.1.45 6
Ke_PR6= 0.24;  KmGOA_PR6=0.15;  KmSER_PR6=2.7;  KmGLY_PR6=33.0;
%1.1.1.29 7
Ke_PR7= 250000.0;  KiHPR_PR7= 12.0;  KmHPR_PR7=0.09;
%2.7.1.31 8
Ke_PR8= 300.0;  KmATP_PR8= 0.21;  KmGCEA_PR8=0.25;  KiPGA_PR8=0.72;%KmATP_PR8= 0.21;  KmGCEA_PR8=0.25;  KiPGA_PR8=0.36;
%Tgca 9
KmGCA_PR9= 0.2;  KiGCEA_PR9= 0.22;
%Tgcea 10
KmGCEA_PR10= 0.39;  KiGCA_PR10= 0.28;

% Transport coeffcient (1/second)
Voaa =1.5;
Vmal =1.5;
Vpyr =1.5;
Vpep =1.5;
Vt =1.5;
Vleak=1;
Vpga=2;

KmPGA_62 = 0.08; 
KmPEP_62 = 0.3;
Ke_62=0.4302;%G66 = +0.5; 



global KValue;
KValue=zeros(48,10);
KValue(1,1)=KmCO2_1;  KValue(1,2)=Ke_1;
KValue(2,1)=KmHCO3_2;  KValue(2,2)= KmPEP_2;   KValue(2,3)=Kimal_2;
KValue(3,1)=KmNADPH_3;  KValue(3,2)=KmOAA_3;  KValue(3,3)=KmNADP_3;  KValue(3,4)=Kmmal_3;  KValue(3,5)=Ke_3;
KValue(4,1)=KmCO2_4;  KValue(4,2)=KmNADP_4;  KValue(4,3)=KmNADPH_4;  KValue(4,4)=KmPyr_4;  KValue(4,5)=Kmmal_4;  KValue(4,6)=Ke_4;
KValue(5,1)=KiPEP_5;  KValue(5,2)=KmATP_5;  KValue(5,3)=KmPyr_5;
KValue(6,1)=KmCO2_6;  KValue(6,2)=KmO2_6;  KValue(6,3)=KmRuBP_6;  KValue(6,4)=KiPGA_6;  KValue(6,5)=KiFBP_6;  KValue(6,6)=KiSBP_6;  KValue(6,7)=KiPi_6;  KValue(6,8)=KiNADPH_6;
KValue(7,1)=KmADP_7;  KValue(7,2)=KmATP_7;  KValue(7,3)=KmPGA_7;
KValue(8,1)=KmDPGA_8;  KValue(8,2)=KmNADPH_8;
KValue(9,1)=Ke_9; 
KValue(10,1)=KmDHAP_10;  KValue(10,2)=KmFBP_10;  KValue(10,3)=KmGAP_10;  KValue(10,4)=Ke_10;
KValue(11,1)=KiF6P_11;  KValue(11,2)=KiPi_11;  KValue(11,3)=KmFBP_11;  KValue(11,4)=Ke_11;
KValue(12,1)=KmDHAP_12;  KValue(12,2)=KmE4P_12;  KValue(12,3)=Ke_12;
KValue(13,1)=KiPi_13;  KValue(13,2)=KmSBP_13;  KValue(13,3)=Ke_13;
KValue(14,1)=KmE4P_14;  KValue(14,2)=KmF6P_14;  KValue(14,3)=KmGAP_14;  KValue(14,4)=KmXu5P;  KValue(14,5)=Ke_14;
KValue(15,1)=KmGAP_15;  KValue(15,2)=KmRi5P_15;  KValue(15,3)=KmS7P_15;  KValue(15,4)=KmXu5P_15;  KValue(15,5)=Ke_15;
KValue(16,1)=Ke_16;
KValue(17,1)=Ke_17;
KValue(18,1)=KiADP_18;  KValue(18,2)=Ki_ADP_18;  KValue(18,3)=KiPGA_18;  KValue(18,4)=KiPi_18;  KValue(18,5)=KiRuBP_18;  KValue(18,6)=KmATP_18;  KValue(18,7)=KmRu5P_18;  KValue(18,8)=Ke_18;

KValue(19,1)=KmADP_7Mchl;  KValue(19,2)=KmATP_7Mchl;  KValue(19,3)=KmPGA_7Mchl;
KValue(20,1)=KmDPGA_8Mchl;  KValue(20,1)=KmNADPH_8Mchl;

KValue(21,1)=KiADP_Starch;  KValue(21,2)=KmATP_Starch;  KValue(21,3)=KmG1P_Starch;  KValue(21,4)=KaF6P_Starch;  KValue(21,5)=KaFBP_Starch;  KValue(21,6)=KaPGA_Starch;  KValue(21,7)=Ke_Starch1;   KValue(21,8)=Ke_Starch2;
KValue(22,1)=KmPGA_PGASink;
KValue(23,1)=KmDHAP_Suc1;  KValue(23,2)=KmGAP_Suc1;  KValue(23,3)=KmFBP_Suc1;  KValue(23,4)=Ke_Suc1;
KValue(24,1)=KiF26BP_Suc2;  KValue(24,2)=KiF6P_Suc2;  KValue(24,3)=KiPi_Suc2;  KValue(24,4)=KmFBP_Suc2;  KValue(24,5)=Ke_Suc2;
KValue(25,1)=Ke_Suc5;  KValue(25,2)=Ke_Suc6;
KValue(26,1)=KmG1P_Suc7;  KValue(26,2)=KmPPi_Suc7;  KValue(26,3)=KmUDPG_Suc7;  KValue(26,4)=KmUTP_Suc7;  KValue(26,5)=Ke_Suc7;
KValue(27,1)=KiFBP_Suc8;  KValue(27,2)=KiPi_Suc8;  KValue(27,3)=KiSuc_Suc8;  KValue(27,4)=KiSucP_Suc8;  KValue(27,5)=KiUDP_Suc8;  KValue(27,6)=KmF6P_Suc8;  KValue(27,7)=KmUDPG_Suc8;  KValue(27,8)=Ke_Suc8; 
KValue(28,1)=KmSuc_Suc9;  KValue(28,2)=KmSucP_Suc9;  KValue(28,3)=Ke_Suc9;
KValue(29,1)=KmSuc_Suc10;
KValue(30,1)=KiADP_Suc3;  KValue(30,2)=KIDHAP_Suc3;  KValue(30,3)=KmATP_Suc3;  KValue(30,4)=KmF26BP_Suc3;  KValue(30,5)=KmF6P_Suc3;  KValue(30,6)=Ke_Suc3;
KValue(31,1)=KiF6P_Suc4;  KValue(31,2)=KiPi_Suc4;  KValue(31,3)=KmF26BP_Suc4;  
KValue(36,1)=KePi;

KValue(32,1)=KmADP_ATPM;  KValue(32,2)=KmATP_ATPM;  KValue(32,3)=KmPi_ATPM;  KValue(32,4)=X;  KValue(32,5)=Y;  KValue(32,6)=F;  KValue(32,7)=Q;  KValue(32,8)=D; KValue(32,9)=Ke_ATPM;
KValue(33,1)=KmNADP_NADPHM;  KValue(33,2)=KmNADPH_NADPHM; KValue(33,3)=Ke_NADPHM;  KValue(33,4)=E; 
KValue(34,1)=KmADP_ATPB;  KValue(34,2)=KmPi_ATPB;   KValue(34,3)=KmATP_ATPB;  KValue(34,4)=Ke_ATPB;   KValue(34,5)=G;
KValue(37,1)=KmNADP_NADPHB;  KValue(37,2)=KmNADPH_NADPHB; KValue(37,3)=Ke_NADPHB; 

KValue(35,1)=Voaa;  KValue(35,2)=Vmal;  KValue(35,3)=Vpyr;  KValue(35,4)=Vpep;  KValue(35,5)=Vt;  KValue(35,6)=Vleak; KValue(35,7)=Vpga;


KValue(38,1)= KmCO2_PR1; KValue(38,2)= KmO2_PR1;  KValue(38,3)=KmRuBP_PR1;  KValue(38,4)=KiPGA_PR1;  KValue(38,5)=KiFBP_PR1;  KValue(38,6)=KiSBP_PR1;  KValue(38,7)=KiPi_PR1;  KValue(38,8)=KiNADPH_PR1;
KValue(39,1)=KmPGCA_PR2;  KValue(39,2)=KiPI_PR2;  KValue(39,3)=KiGCA_PR2;
KValue(40,1)=KmGCA_PR3;
KValue(41,1)=Ke_PS4;  KValue(41,2)=KmGOA_PS4;  KValue(41,3)=KmGLU_PS4;  KValue(41,4)=KiGLY_PS4;
KValue(42,1)=KmGLY_PS5;  KValue(42,2)=KiSER_PS5;
KValue(43,1)=Ke_PR6;  KValue(43,2)=KmGOA_PR6;  KValue(43,3)=KmSER_PR6;  KValue(43,4)=KmGLY_PR6;
KValue(44,1)=Ke_PR7;  KValue(44,2)=KiHPR_PR7;  KValue(44,3)=KmHPR_PR7;
KValue(45,1)=Ke_PR8;  KValue(45,2)=KmATP_PR8;  KValue(45,3)=KmGCEA_PR8;  KValue(45,4)=KiPGA_PR8;
KValue(46,1)=KmGCA_PR9;  KValue(46,2)=KiGCEA_PR9;
KValue(47,1)=KmGCEA_PR10;  KValue(47,2)=KiGCA_PR10;

KValue(48,1)=KmPGA_62; KValue(48,2)=KmPEP_62; KValue(48,3)=Ke_62;

%mM/(L*s) 
global Enzyme;
global Vpmax;
global Vcmax;
global MeasuredTemp;
global FactorVP;
global FactorVC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Temp correction
Ea_Vpmax=94.8*1000;
dS_Vpmax=0.25*1000;
Hd_Vpmax=73.3*1000;
Ea_Vcmax=78*1000;
Ea_PPDK=58.1*1000;
R=8.3144598;% m2 kg s-2 K-1 mol-1
MTempCorr_V2=exp(Ea_Vpmax *((MeasuredTemp+273.15)-298.15)/(298.15*R*(MeasuredTemp+273.15)))*(1+exp((298.15*dS_Vpmax-Hd_Vpmax)/(298.15*R)))/(1+exp(((MeasuredTemp+273.15)*dS_Vpmax-Hd_Vpmax)/((MeasuredTemp+273.15)*R)));
MTempCorr_V6=exp(Ea_Vcmax*((MeasuredTemp+273.15)-298.15)/(298.15*R*(MeasuredTemp+273.15)));
MTempCorr_V5=exp(Ea_PPDK*((MeasuredTemp+273.15)-298.15)/(298.15*R*(MeasuredTemp+273.15)));
Vpmax_25=Vpmax/MTempCorr_V2;
Vcmax_25=Vcmax/MTempCorr_V6;
Vppdkmax_25=Vcmax/MTempCorr_V5;
Vm_2 = Vpmax_25/1000/FactorVP;%3.6;%4; 
Vm_1 = 200;
Vm_6o = 0.065;%2.4;%1.8;%2.9;%3.4967170968;2.913930914 2.09803025808
Vppdkmax_E=Vppdkmax_25/1000;
Vm_6 = Vcmax_25/1000/FactorVC;
Vm_3 = 1.8*Vppdkmax_E*VMDHsen;%/2.4*Vm_2;%/1.2;        
Vm_4 = 1.33*Vppdkmax_E*VMEsen;%/2.4*Vm_2;%/1.5;%2.4;%3.2
Vm_5 = 1.33*Vppdkmax_E*VPPDKsen;%/2.4*Vm_2;%*1.5;%2.4; WY 201902
Vm_78= 0.4*Vm_6/Vm_6o*VGAPDHsen;%Vm_7 = 30.15;%16.07510272;% 30.15;%16.07510272
Vm_8 = 0;%4.308781695;
Vm_10 = 0.0731*Vm_6/Vm_6o;%2.9253370968*1.5;%1.4626685484;
Vm_11 =0.0436*Vm_6/Vm_6o*VFBPsen;
Vm_12 = 0.1097*Vm_6/Vm_6o;
Vm_13 = 0.0292*Vm_6/Vm_6o*VSBPsen;
Vm_14 = 0.2810*Vm_6/Vm_6o;
Vm_15 = 0.2810*Vm_6/Vm_6o;
Vm_18 = 1.7552*Vm_6/Vm_6o*VPRKsen;


Vm_78Mchl=0.3*Vm_6/Vm_6o*VGAPDHsen;%Vm_7Mchl = 15.1;
Vm_8Mchl =0;%2.6929885593*Vm_6/Vm_6o;
Vm_Starch =0.0133*Vm_6/Vm_6o;
Vm_Sta1=0.03*Vm_6/Vm_6o;
Vm_Sta2=1*Vm_6/Vm_6o;
Vm_Sta3=0.025*Vm_6/Vm_6o;

Vm_PGASink = 0.002*Vm_6/Vm_6o;%0.5/5;
Vm_Suc1 = 0.0081*Vm_6/Vm_6o;
Vm_Suc2 = 0.0064*Vm_6/Vm_6o;
Vm_Suc7 = 0.0058*Vm_6/Vm_6o;
Vm_Suc8 = 0.0278*Vm_6/Vm_6o;
Vm_Suc9 = 0.0278*Vm_6/Vm_6o;
Vm_Suc10 =0.0035*Vm_6/Vm_6o; %2.0
Vm_Suc3 =0.0010*Vm_6/Vm_6o;
Vm_Suc4 =8.4096e-04*Vm_6/Vm_6o;

Jmax =0.5*Vm_6/Vm_6o*Jmaxsen;% 20;
Vm_ATPM = 0.3*Vm_6/Vm_6o;
Vm_NADPHM = 0.2*Vm_6/Vm_6o; 
Vm_ATPB = 0.3*Vm_6/Vm_6o;
Vm_NADPHB= 0.2*Vm_6/Vm_6o;

Vm_PR1= Vm_6*0.11;%0.69934341936;(Cousins 2010 0.11)
Vm_PR2= 2.6210*Vm_6/Vm_6o;
Vm_PR3= 0.0728*Vm_6/Vm_6o;
Vm_PR4= 0.1373*Vm_6/Vm_6o;
Vm_PR5= 0.1247*Vm_6/Vm_6o;
Vm_PR6= 0.1653*Vm_6/Vm_6o;
Vm_PR7= 0.5005*Vm_6/Vm_6o;
Vm_PR8= 0.2858*Vm_6/Vm_6o;
VTgca_PR9=0.3*Vm_6/Vm_6o;
VTgcea_PR10=0.25*Vm_6/Vm_6o;

Vm_62 =0.001*Vexsen;%mM/s %1.45 E-5; 

Vm_OAA_M=0.08*Vm_6/Vm_6o;
Vm_PYR_B=0.15*Vm_6/Vm_6o;
Vm_PYR_M=0.15*Vm_6/Vm_6o;
Vm_PEP_M=0.15*Vm_6/Vm_6o;
Vtp_Bchl=0.75;
Vtp_Mchl=0.75;

% transport between two cell types
global phi;
global Lpd;

Pmal=0.0421*(phi/0.03)/(Lpd/400);
Ppyr=0.0436*(phi/0.03)/(Lpd/400);
Pco2=0.1139*(phi/0.03)/(Lpd/400);
PC3P=0.0327;
Pc3p=PC3P*(phi/0.03)/(Lpd/400);
Pco2_B=0.4;%%PCO2_B=0.002 cm s-1 SChl/Sl =10

global MRd;
global vrpd;
vrpd=MRd/2/1000;%0.0005;

global Velocity_s;
Velocity_s=zeros(58,1); 
Velocity_s(1)=Vm_1;
Velocity_s(2)=Vm_2;          
Velocity_s(3)=Vm_3;
Velocity_s(4)=Vm_4;        
Velocity_s(5)=Vm_5;
Velocity_s(6)=Vm_6;
Velocity_s(7)=Vm_78;
Velocity_s(8)=Vm_8;
Velocity_s(9)=Vm_10;
Velocity_s(10)=Vm_11;
Velocity_s(11)=Vm_12;
Velocity_s(12)=Vm_13;
Velocity_s(13)=Vm_14;
Velocity_s(14)=Vm_15;
Velocity_s(15)=Vm_18;
Velocity_s(16)=Vm_78Mchl;
Velocity_s(17)=Vm_8Mchl;
Velocity_s(18)=Vm_Starch;
Velocity_s(19)=Vm_PGASink;
Velocity_s(20)=Vm_Suc1;
Velocity_s(21)=Vm_Suc2;
Velocity_s(22)=Vm_Suc7;
Velocity_s(23)=Vm_Suc8;
Velocity_s(24)=Vm_Suc9;
Velocity_s(25)=Vm_Suc10;
Velocity_s(26)=Vm_Suc3;
Velocity_s(27)=Vm_Suc4;
Velocity_s(28)=0;%Radiation_PAR*Convert/1000;
Velocity_s(29)=Jmax;
Velocity_s(30)=Vm_ATPM;
Velocity_s(31)=Vm_NADPHM; 
Velocity_s(32)=Vm_ATPB;
Velocity_s(33)=Vm_NADPHB;
Velocity_s(34)=Vm_PR1;
Velocity_s(35)=Vm_PR2;
Velocity_s(36)=Vm_PR3;
Velocity_s(37)=Vm_PR4;
Velocity_s(38)=Vm_PR5;
Velocity_s(39)=Vm_PR6;
Velocity_s(40)=Vm_PR7;
Velocity_s(41)=Vm_PR8;
Velocity_s(42)=VTgca_PR9;
Velocity_s(43)=VTgcea_PR10;
Velocity_s(44)=Vm_62;
Velocity_s(45)=Vtp_Bchl;
Velocity_s(46)=Vtp_Mchl;
Velocity_s(47)=Vm_Sta1;
Velocity_s(48)=Vm_Sta2;
Velocity_s(49)=Vm_Sta3;
Velocity_s(50)=Vm_OAA_M;
Velocity_s(51)=Vm_PYR_B;
Velocity_s(52)=Vm_PYR_M;
Velocity_s(53)=Vm_PEP_M;
Velocity_s(54)=Pmal;
Velocity_s(55)=Ppyr;
Velocity_s(56)=Pco2;
Velocity_s(57)=Pc3p;
Velocity_s(58)=Pco2_B;

 
