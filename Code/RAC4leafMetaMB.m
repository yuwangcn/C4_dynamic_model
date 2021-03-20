function LM_DYDT = RAC4leafMetaMB(t,LM_con)  
%global Para_mata;
global Jsen;
global Radiation_PARo;
global Radiation_PAR;
global vrpd; 
global WeatherTemperature;
Convert=1E6/(2.35E5); 
t00=300;
if t<t00
    Radiation_PAR=00/Convert*Jsen;
end
if t>=t00
Radiation_PAR=Radiation_PARo/Convert*Jsen;
end
% if t>=t00+1800&&t<t00+3600
%     Radiation_PAR=200/Convert*Jsen;
% end
% if t>=t00+3600&&t<t00+5400
%     Radiation_PAR=1800/Convert*Jsen;
% end

%RuACT_Con=LM_con(95:98);
global PSPR_RA_O2;
global PSPR_RA_CO2;
PSPR_RA_CO2=LM_con(37);
PSPR_RA_O2=LM_con(81);
%global PSPR_RA_CA;
global PS2RA_ATP;
PS2RA_ATP=LM_con(41);
RuACT_mb	=	zeros(4,1)	;				

LM_v=RAC4LeafMetaVel(t,LM_con);
NetAssimilation=LM_v(1);
vCO2b=LM_v(2);
vCO2s=LM_v(3);
vH2Ob=LM_v(4);
vH2Os=LM_v(5);
EnergyBalanceResidual=LM_v(6);
vCO2total=LM_v(7);
vH2Ototal=LM_v(8);
vgs=LM_v(9);
Enz_v=LM_v(10:111,1);
v1=Enz_v(1);
v2=Enz_v(2);
v3=Enz_v(3);
v4=Enz_v(4);
v5=Enz_v(5);
v6=Enz_v(6);
v7=Enz_v(7);
v8=Enz_v(8);
v10=Enz_v(9);
v11=Enz_v(10);
v12=Enz_v(11);
v13=Enz_v(12);
v14=Enz_v(13);
v15=Enz_v(14);
v18=Enz_v(15);
v7Mchl=Enz_v(16);
v8Mchl=Enz_v(17);
vStarch1=Enz_v(18);
vPGASink=Enz_v(19);
vSuc1=Enz_v(20);
vSuc2=Enz_v(21);
vSuc3=Enz_v(22);
vSuc4=Enz_v(23);
vSuc7=Enz_v(24);
vSuc8=Enz_v(25);
vSuc9=Enz_v(26);
vSuc10=Enz_v(27);
vATPM=Enz_v(28);
vNADPHM=Enz_v(29);
vATPB=Enz_v(30);
vOAA_M=Enz_v(31);
vMAL_M=Enz_v(32);
vMAL=Enz_v(33);
vMAL_B=Enz_v(34);
vPYR_B=Enz_v(35);
vPYR=Enz_v(36);
vPYR_M=Enz_v(37);
vPEP_M=Enz_v(38);
vPGA_B=Enz_v(39);
vPGA=Enz_v(40);
vPGA_M=Enz_v(41);
vGAP_M=Enz_v(42);
vGAP=Enz_v(43);
vGAP_B=Enz_v(44);
vDHAP_M=Enz_v(45);
vDHAP=Enz_v(46);
vDHAP_B=Enz_v(47);
vleak_B=Enz_v(48);
vleak=Enz_v(49);
vNADPHB=Enz_v(50);
vpr1=Enz_v(51);
vpr2=Enz_v(52);
vpr3=Enz_v(53);
vpr4=Enz_v(54);
vpr5=Enz_v(55);
vpr6=Enz_v(56);
vpr7=Enz_v(57);
vpr8=Enz_v(58);
vpr9=Enz_v(59);
vpr10=Enz_v(60);
vinf=Enz_v(61);
v62=Enz_v(62);
v78=Enz_v(63);
v78Mchl=Enz_v(64);
vStarch2=Enz_v(65);
vgly1B=Enz_v(66);
vhexp=Enz_v(67);
vSta1=Enz_v(68);
vSta2=Enz_v(69);
vSta3=Enz_v(70);
vtATP=Enz_v(71);

vPCK1=Enz_v(72);
vPCK2=Enz_v(73);
vPCK3=Enz_v(74);
vPCK4=Enz_v(75);
vPCK5=Enz_v(76);
vPCK6=Enz_v(77);
vAsp=Enz_v(78);
vAla=Enz_v(79);
vPEP=Enz_v(80);
vOAA_B=Enz_v(81);
vATP_B=Enz_v(82);

vO2_Mchl=Enz_v(83);
vO2_Bchl=Enz_v(84);
vtO2=Enz_v(85);
vtO2_B=Enz_v(86);
vtO2_M=Enz_v(87);
v5B=Enz_v(88);
vPEP_B=Enz_v(89);
vpr10M=Enz_v(90);
vpr8M=Enz_v(91);
vMDH_BM=Enz_v(92);
vNADME=Enz_v(93);
vtOAA_Bm=Enz_v(94);
vPYR_Bm=Enz_v(95);
vMAL_Bm=Enz_v(96);
vNAE=0;
vAspB=Enz_v(97);
vAlaB=Enz_v(98);
vAspM=Enz_v(99);
vAlaM=Enz_v(100);
vEA_PPDKRP_I=Enz_v(101);
vEA_PPDKRP_A=Enz_v(102);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EnzAct_v=LM_v(112:123,1);
vATPsynthase_Act_Mchl=EnzAct_v(1);
vNADPMDH_Act=EnzAct_v(2);
vGAPDH_Act_Mchl=EnzAct_v(3);
vATPsynthase_Act_Bchl=EnzAct_v(4);
vPEPC_Act=EnzAct_v(5);
ActRubisco0=EnzAct_v(6);
vRubisco_Act=EnzAct_v(7);
vGAPDH_Act_Bchl=EnzAct_v(8);
vFBPase_Act=EnzAct_v(9);
vSBPase_Act=EnzAct_v(10);
vPRK_Act=EnzAct_v(11);
vRCA_Act=EnzAct_v(12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VolMC=0.01;
VolMchl=0.02;
VolBSC=0.0045;
VolBchl=0.009;
Volper=0.00045;
Volmito=0.00045;

global option2; %%O2 diffusion
option2=1;
   
%Delta_MC_CO2=(vinf-v1+vleak+2*vrpd)/VolMC;
Delta_MC_CO2=(vinf-v1+vleak+vrpd)/VolMC;
Delta_MC_HCO3= (v1 - v2)/VolMC;
Delta_MC_Malate= (vMAL_M - vMAL)/VolMC;
Delta_MC_PGA= (vPGA - vPGA_M -v62- vPGASink)/VolMC;
Delta_MC_FBP= (vSuc1- vSuc2)/VolMC;
Delta_MC_UDPG= (vSuc7 - vSuc8)/VolMC;
Delta_MC_SUCP= (vSuc8 - vSuc9)/VolMC;
Delta_MC_SUC= (vSuc9 - vSuc10)/VolMC;
Delta_MC_F26BP= (vSuc3 - vSuc4)/VolMC;
Delta_MC_ATP= (vtATP-vSuc7-vSuc3)/VolMC;
Delta_MC_T3P= (vGAP_M + vDHAP_M - vGAP - vDHAP - 2*vSuc1)/VolMC;
Delta_MC_HexP= (vSuc2 + vSuc4 - vSuc3 - vSuc7 - vSuc8)/VolMC;
Delta_MC_Sucrose= (vSuc10)/VolMC;

Delta_Mchl_OAA= (vOAA_M - v3)/VolMchl;
Delta_Mchl_Malate= (v3 - vMAL_M)/VolMchl;
Delta_Mchl_PEP= (v5 - vPEP_M)/VolMchl;
Delta_Mchl_Pyruvate= (vPYR_M - v5)/VolMchl;
Delta_Mchl_NADPH= (vNADPHM - v3 - v78Mchl)/VolMchl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mchl_TotalATPsynthase=0.00545;
Mchl_TotalGAPDH=0.00492;
Mchl_TotalNADPMDH=0.00006;

if vATPsynthase_Act_Mchl>=0
    vATPsynthase_Act_MchlF=vATPsynthase_Act_Mchl;
else
    vATPsynthase_Act_MchlF=0;
end
if vNADPMDH_Act>=0
    vNADPMDH_ActF=vNADPMDH_Act;
else
    vNADPMDH_ActF=0;
end
if vGAPDH_Act_Mchl>=0
    vGAPDH_Act_MchlF=vGAPDH_Act_Mchl;
else
    vGAPDH_Act_MchlF=0;
end
Delta_Mchl_NADPH= (vNADPHM - v3 - v78Mchl-vATPsynthase_Act_MchlF*Mchl_TotalATPsynthase-vNADPMDH_ActF*Mchl_TotalNADPMDH-vGAPDH_Act_MchlF*Mchl_TotalGAPDH)/VolMchl;%WY1911
Delta_Mchl_ATP= (vATPM - 2*v5 - v78Mchl-vtATP - vpr8M)/VolMchl;
Delta_Mchl_PGA= (vPGA_M - v78Mchl)/VolMchl;
Delta_Mchl_DPGA= 0;
Delta_Mchl_T3P= (v78Mchl - vGAP_M - vDHAP_M)/VolMchl;

Delta_BSC_T3P= (vGAP + vDHAP - vGAP_B - vDHAP_B)/VolBSC;
Delta_BSC_PGA= (vPGA_B - vPGA)/VolBSC;

Delta_Bchl_CO2= (v4 - v6 - vleak_B)/VolBchl;
Delta_Bchl_RuBP= (v18 - v6-vpr1)/VolBchl;
Delta_Bchl_PGA= (2*v6 - v78 - vPGA_B + vpr1 + vpr8)/VolBchl;
Delta_Bchl_DPGA=0;
Delta_Bchl_SBP= (v12 - v13)/VolBchl;
Delta_Bchl_S7P= (v13 - v15)/VolBchl;
Delta_Bchl_FBP= (v10 - v11+vgly1B)/VolBchl;
Delta_Bchl_E4P= (v14 - v12)/VolBchl;
Delta_Bchl_Starch= vSta3/VolBchl;
Delta_Bchl_Rubisco= 0;
Delta_Bchl_T3P= (vGAP_B + vDHAP_B + v78 - 2*v10 - v14 - v15 - v12)/VolBchl;
Delta_Bchl_HexP= (v11 - v14 - vSta1 -vgly1B+vhexp)/VolBchl;
Delta_Bchl_Pent= (v14 + 2*v15 - v18)/VolBchl;
Delta_Bchl_Pyruvate= (v4 - vPYR_B-v5B)/VolBchl;%%%%WY
Delta_Bchl_PPi= (vSta1-vSta2)/VolBchl;
Delta_Bchl_ADPG= (vSta1-vSta3)/VolBchl;
Delta_Bchl_PGCA= (vpr1 - vpr2)/VolBchl;
Delta_Bchl_GCA= (vpr2 - vpr9)/VolBchl;
Delta_Bchl_GCEA= (vpr10 - vpr8)/VolBchl;

Delta_Bper_GCA= (vpr9 - vpr3)/Volper;
Delta_Bper_GOA= (vpr3 - vpr4 -vpr6)/Volper;
Delta_Bper_GLY= (vpr4 + vpr6 - 2*vpr5)/Volper;
Delta_Bper_SER= (vpr5 - vpr6)/Volper;
Delta_Bper_HPR= (vpr6 - vpr7)/Volper;
Delta_Bper_GCEA= (vpr7 - vpr10)/Volper;

Delta_MC_O2=0;
Delta_Mchl_O2=vO2_Mchl-vtO2_M;
Delta_BSC_O2=vtO2_B-vtO2;
Delta_Bchl_O2=vO2_Bchl-vtO2_B;

   Delta_MC_OAA= (v2 - vOAA_M-vPCK1)/VolMC;
   Delta_MC_Pyruvate= (vPYR - vPYR_M+vPCK5)/VolMC;
   Delta_MC_PEP= (vPEP_M - v2 + v62 + vPEP )/VolMC;

   Delta_BSC_PEP=(vPCK3-vPEP+vPEP_B)/VolBSC;%%%%%%WY

   if vATPsynthase_Act_Bchl>=0
    vATPsynthase_Act_BchlF=vATPsynthase_Act_Bchl;
    else
    vATPsynthase_Act_BchlF=0;
    end
%     if vPEPC_Act>=0
%     vPEPC_ActF=vPEPC_Act;
%     else
%     vPEPC_ActF=0;
%     end
    if vGAPDH_Act_Bchl>=0
    vGAPDH_Act_BchlF=vGAPDH_Act_Bchl;
    else
    vGAPDH_Act_BchlF=0;
    end
    if vFBPase_Act>=0
    vFBPase_ActF=vFBPase_Act;
    else
    vFBPase_ActF=0;
    end
    if vSBPase_Act>=0
    vSBPase_ActF=vSBPase_Act;
    else
    vSBPase_ActF=0;
    end
    if vPRK_Act>=0
    vPRK_ActF=vPRK_Act;
    else
    vPRK_ActF=0;
    end
    Bchl_TotalATPsynthase=0.00545;
    Bchl_TotalRca=0.00160;
    Bchl_TotalGAPDH=0.00328;
    Bchl_TotalFBPase=0.00286;
    Bchl_TotalSBPase=0.0006;
    Bchl_TotalPRK=0.00130;

   Delta_Bchl_NADPH=(v4 - v78-vPCK6+vNADPHB-vATPsynthase_Act_BchlF*Bchl_TotalATPsynthase-vGAPDH_Act_BchlF*Bchl_TotalGAPDH-vFBPase_ActF*Bchl_TotalFBPase-vSBPase_ActF*Bchl_TotalSBPase-vPRK_ActF*Bchl_TotalPRK)/VolBchl;%+ vNADPHB;

   Delta_Bchl_Malate= (vMAL_B - v4+vPCK6)/VolBchl; 
   Delta_Bchl_OAA=(vOAA_B-vPCK6)/VolBchl;
%Delta_Bchl_ATP=(vATPB - v78 - v18 - vSta1 - vpr8-vgly1B-vATP_B-2*v5B)/VolBchl;%%%%%%WY
    if vRubisco_Act>=0
    vRubisco_ActF=vRubisco_Act;
    else
    vRubisco_ActF=0;
    end
   Bchl_TotalRubisco=0.0159;
   Delta_Bchl_ATP=(vATPB - v78 - v18 - vSta1 - vpr8-vgly1B-vATP_B-2*v5B-vRubisco_ActF*Bchl_TotalRubisco)/VolBchl;%%%%%%WY1911
   Delta_MC_Glu=(vPCK5-vPCK1)/VolMC;
   Delta_MC_OxoG=(vPCK1-vPCK5)/VolMC;
   Delta_MC_Asp=(vPCK1-vAsp)/VolMC;
   Delta_MC_Ala=(vAla-vPCK5)/VolMC;
   Delta_BSC_OxoG=(vPCK4-vPCK2)/VolBSC;
   Delta_BSC_Glu=(vPCK2-vPCK4)/VolBSC;
   Delta_BSC_Asp=(vAsp-vPCK2)/VolBSC;
   Delta_BSC_Ala=(vPCK4-vAla)/VolBSC;
   Delta_Bchl_PEP=(v5B-vPEP_B)/VolBchl;
   Delta_Mchl_GCEA=(vpr10M - vpr8M)/VolMchl;  
   
%    Delta_BSC_CO2= (vpr5+vleak_B - vleak+vrpd+vPCK3)/VolBSC;
%    Delta_BSC_ATP=vATP_B-vPCK3;
%    Delta_BSC_Pyruvate= (vPYR_B - vPYR-vPCK4)/VolBSC;
%    Delta_BSC_OAA=(vPCK2-vOAA_B-vPCK3)/VolBSC;
%    Delta_BSC_Malate= (vMAL - vMAL_B)/VolBSC;

   
   
 Delta_BSC_CO2= (vpr5+vleak_B - vleak+vrpd+vPCK3+vNADME)/VolBSC;
%   Delta_BSC_CO2= (vpr5+vleak_B - vleak+vPCK3+vNADME)/VolBSC;
   Delta_BSC_ATP=(vATP_B-vPCK3+vNAE*2.5)/VolBSC;
   Delta_BSC_Pyruvate= (vPYR_B - vPYR-vPCK4+vPYR_Bm)/VolBSC;
   Delta_BSC_OAA=(vPCK2-vOAA_B-vPCK3-vtOAA_Bm)/VolBSC;
   Delta_BSC_Malate= (vMAL - vMAL_B-vMAL_Bm)/VolBSC;
   
   Delta_Bmito_OAA=(vtOAA_Bm-vMDH_BM)/Volmito;
   Delta_Bmito_MAL=(vMDH_BM+vMAL_Bm-vNADME)/Volmito;
   Delta_Bmito_PYR=(vNADME-vPYR_Bm)/Volmito;
   Delta_Bmito_CO2=0;
   Delta_Bmito_NADH=(vNADME-vMDH_BM-vNAE)/Volmito;
      Delta_Bchl_Asp=vAspB;
   Delta_Bchl_Ala=vAlaB;
   Delta_Mchl_Asp=vAspM;
   Delta_Mchl_Ala=vAlaM;
   
 Delta_E_PPDKdt=(vEA_PPDKRP_A-vEA_PPDKRP_I)/VolMchl;
 Delta_EP_PPDKdt=(vEA_PPDKRP_I-vEA_PPDKRP_A)/VolMchl;
if option2 == 0
    Delta_MC_O2=0;
    Delta_Mchl_O2=0;
    Delta_BSC_O2=0;
    Delta_Bchl_O2=0;
end


Cpwater=4.184;%Jg-1c-1
Vandp=198;%gm-2 leaf thickness 200um leaf density 0.7*10^3kg m-3
%leaf volume leaf thickness 200um~0.2L
%bundary layer
%leaf air space 20%
Vol_airspace=0.04;%L
Molar_Volume=22.4/273*(WeatherTemperature+273);
Delta_Ci=(vCO2s-NetAssimilation)/(Vol_airspace/Molar_Volume);
Delta_Ci=(vCO2total-NetAssimilation)/(Vol_airspace/Molar_Volume);% WY202010

Delta_Cb=0;%vCO2b-vCO2s;
Delta_Eb=0;%vH2Os-vH2Ob;
Delta_Gs=vgs;
Delta_Tleaf=EnergyBalanceResidual/Cpwater/Vandp;
Delta_H2Oou=vH2Ob/10^6;%umol->mol
Delta_CO2in=vCO2b/10^6;
Delta_H2Oou=vH2Ototal/10^6;%umol->mol %WY202010
Delta_CO2in=vCO2total/10^6;

LeafMB=zeros(7,1);
LeafMB(1)=Delta_Ci;
LeafMB(2)=Delta_Cb;
LeafMB(3)=Delta_Eb;
LeafMB(4)=Delta_Gs;
LeafMB(5)=Delta_Tleaf;
LeafMB(6)=Delta_H2Oou;
LeafMB(7)=Delta_CO2in;

% %%%%%%%%%%%%%%%%%%%%%%%%
% %Redox regulaltion ODEs
% %Delta_Mchl_ActThioredoxin=(vReTRX_M-vReATPsyn_M-vReGAPDH_M-vReNADPMDH)/VolMchl;
% Delta_Mchl_ActATPsynthase=(vReATPsyn_M)/VolMchl;
% Delta_Mchl_ActGAPDH=(vReGAPDH_M)/VolMchl;
% Delta_Mchl_ActNADPMDH=(vReNADPMDH)/VolMchl;
% 
% %Delta_Bchl_ActThioredoxin=(vReTRX_B-vReATPsyn_B-vReRca_B-vReGAPDH_B+vReFBPase_B+vReSBPase_B+vRePRK_B)/VolBchl;
% Delta_Bchl_ActATPsynthase=(vReATPsyn_B)/VolBchl;
% Delta_Bchl_ActRca=(vReRca_B)/VolBchl;
% Delta_Bchl_ActGAPDH=(vReGAPDH_B)/VolBchl;
% Delta_Bchl_ActFBPase=(vReFBPase_B)/VolBchl;
% Delta_Bchl_ActSBPase=(vReSBPase_B)/VolBchl;
% Delta_Bchl_ActPRK=(vRePRK_B)/VolBchl;
% Delta_Bchl_ActRubisco=(vPiRubisco)/VolBchl;
% 
% %Enzyme activation 
%%%%%%%%%%%%%%%%%%%%%%%
AEMB=zeros(10,1);

Delta_Mchl_ActATPsynthase=vATPsynthase_Act_Mchl;
Delta_Mchl_ActGAPDH=vGAPDH_Act_Mchl;
Delta_Mchl_ActNADPMDH=vNADPMDH_Act;
Delta_Bchl_ActATPsynthase=vATPsynthase_Act_Bchl;
Delta_Bchl_ActPEPC=vPEPC_Act;
Delta_Bchl_ActGAPDH=vGAPDH_Act_Bchl;
Delta_Bchl_ActFBPase=vFBPase_Act;
Delta_Bchl_ActSBPase=vSBPase_Act;
Delta_Bchl_ActPRK=vPRK_Act;
Delta_Bchl_ActRubisco=vRubisco_Act;
Delta_Bchl_ActRCA=vRCA_Act;

AEMB(1)=Delta_Mchl_ActATPsynthase;
AEMB(2)=Delta_Mchl_ActGAPDH;
AEMB(3)=Delta_Mchl_ActNADPMDH;
AEMB(4)=Delta_Bchl_ActATPsynthase;
AEMB(5)=Delta_Bchl_ActPEPC;
AEMB(6)=Delta_Bchl_ActGAPDH;
AEMB(7)=Delta_Bchl_ActFBPase;
AEMB(8)=Delta_Bchl_ActSBPase;
AEMB(9)=Delta_Bchl_ActPRK;
AEMB(10)=Delta_Bchl_ActRubisco;
AEMB(11)=Delta_Bchl_ActRCA;

EnzMB=zeros(81,1);
EnzMB(1)=Delta_MC_HCO3;
EnzMB(2)=Delta_MC_OAA;
EnzMB(3)=Delta_MC_PEP;
EnzMB(4)=Delta_MC_Malate;
EnzMB(5)=Delta_MC_Pyruvate;
EnzMB(6)=Delta_MC_PGA;
EnzMB(7)=Delta_MC_FBP;
EnzMB(8)=Delta_MC_UDPG;
EnzMB(9)=Delta_MC_SUCP;
EnzMB(10)=Delta_MC_SUC;
EnzMB(11)=Delta_MC_F26BP;
EnzMB(12)=Delta_MC_ATP;
EnzMB(13)=Delta_MC_T3P;
EnzMB(14)=Delta_MC_HexP;
EnzMB(15)=Delta_MC_Sucrose;
EnzMB(16)=Delta_Mchl_OAA;
EnzMB(17)=Delta_Mchl_Malate;
EnzMB(18)=Delta_Mchl_PEP;
EnzMB(19)=Delta_Mchl_Pyruvate;
EnzMB(20)=Delta_Mchl_NADPH;
EnzMB(21)=Delta_Mchl_ATP;
EnzMB(22)=Delta_Mchl_PGA;
EnzMB(23)=Delta_Mchl_DPGA;
EnzMB(24)=Delta_Mchl_T3P;
EnzMB(25)=Delta_BSC_T3P;
EnzMB(26)=Delta_BSC_PGA;
EnzMB(27)=Delta_BSC_Malate;
EnzMB(28)=Delta_BSC_Pyruvate;
EnzMB(29)=Delta_BSC_CO2;
EnzMB(30)=Delta_Bchl_CO2;
EnzMB(31)=Delta_Bchl_RuBP;
EnzMB(32)=Delta_Bchl_PGA;
EnzMB(33)=Delta_Bchl_DPGA;
EnzMB(34)=Delta_Bchl_ATP;
EnzMB(35)=Delta_Bchl_NADPH; 
EnzMB(36)=Delta_Bchl_SBP;
EnzMB(37)=Delta_Bchl_S7P;
EnzMB(38)=Delta_Bchl_FBP;
EnzMB(39)=Delta_Bchl_E4P;
EnzMB(40)=Delta_Bchl_Starch;
EnzMB(41)=Delta_Bchl_Rubisco;
EnzMB(42)=Delta_Bchl_T3P;
EnzMB(43)=Delta_Bchl_HexP;
EnzMB(44)=Delta_Bchl_Pent;
EnzMB(45)=Delta_Bchl_Malate;
EnzMB(46)=Delta_Bchl_Pyruvate;
EnzMB(47)=Delta_Bchl_PGCA;
EnzMB(48)=Delta_Bchl_GCA;
EnzMB(49)=Delta_Bchl_GCEA;
EnzMB(50)=Delta_Bper_GCA;
EnzMB(51)=Delta_Bper_GOA;
EnzMB(52)=Delta_Bper_GLY;
EnzMB(53)=Delta_Bper_SER;
EnzMB(54)=Delta_Bper_HPR;
EnzMB(55)=Delta_Bper_GCEA;
EnzMB(56)=Delta_MC_CO2;
EnzMB(57)=Delta_Bchl_PPi;
EnzMB(58)=Delta_Bchl_ADPG;

EnzMB(59)=Delta_MC_Glu;
EnzMB(60)=Delta_MC_OxoG;
EnzMB(61)=Delta_MC_Asp;
EnzMB(62)=Delta_MC_Ala;

EnzMB(63)=Delta_BSC_OxoG;
EnzMB(64)=Delta_BSC_Glu;
EnzMB(65)=Delta_BSC_Asp;
EnzMB(66)=Delta_BSC_Ala;
EnzMB(67)=Delta_BSC_OAA;
EnzMB(68)=Delta_BSC_PEP;
EnzMB(69)=Delta_BSC_ATP;
EnzMB(70)=Delta_Bchl_OAA;
EnzMB(71)=Delta_MC_O2;
EnzMB(72)=Delta_Mchl_O2;
EnzMB(73)=Delta_BSC_O2;
EnzMB(74)=Delta_Bchl_O2;
EnzMB(75)=Delta_Bchl_PEP;%%%WY PPDK in BS
EnzMB(76)=Delta_Mchl_GCEA;
EnzMB(77)=Delta_Bmito_OAA;
EnzMB(78)=Delta_Bmito_MAL;
EnzMB(79)=Delta_Bmito_PYR;
EnzMB(80)=Delta_Bmito_CO2;
EnzMB(81)=Delta_Bmito_NADH;
EnzMB(82)=Delta_Bchl_Asp;
EnzMB(83)=Delta_Bchl_Ala;
EnzMB(84)=Delta_Mchl_Asp;
EnzMB(85)=Delta_Mchl_Ala;
EnzMB(86)=Delta_E_PPDKdt;
EnzMB(87)=Delta_EP_PPDKdt;

LM_DYDT(1:7,1)=LeafMB;
LM_DYDT(8:94,1)=EnzMB;
LM_DYDT(95:98,1)=RuACT_mb;
LM_DYDT(99:109,1)=AEMB;
%##LM_DYDT(99:109,1)=AEMB;
% %%%%%%%%%%%%%%%%%%%%%%%%
% %Redox regulaltion ODEs
% Delta_Mchl_ActThioredoxin=vReTRX_M-vReATPsyn_M-vReGAPDH_M-vReNADPMDH;
% Delta_Mchl_ActATPsynthase=vReATPsyn_M;
% Delta_Mchl_ActGAPDH=vReGAPDH_M;
% Delta_Mchl_ActNADPMDH=vReNADPMDH;
% Delta_Bchl_ActATPsynthase=vReATPsyn_M;
% Delta_Bchl_ActRca=vReRca_B;
% Delta_Bchl_ActGAPDH=vReGAPDH_B;
% Delta_Bchl_ActFBPase=vReFBPase_B;
% Delta_Bchl_ActSBPase=vReSBPase_B;
% Delta_Bchl_ActPRK=vRePRK_B;
% Delta_Bchl_ActRubisco=vPiRubisco;
% RDMB=zeros(11,1);
% RDMB(1)=Delta_Mchl_ActThioredoxin;
% RDMB(2)=Delta_Mchl_ActATPsynthase;
% RDMB(3)=Delta_Mchl_ActGAPDH;
% RDMB(4)=Delta_Mchl_ActNADPMDH;
% RDMB(5)=Delta_Bchl_ActATPsynthase;
% RDMB(6)=Delta_Bchl_ActRca;
% RDMB(7)=Delta_Bchl_ActGAPDH;
% RDMB(8)=Delta_Bchl_ActFBPase;
% RDMB(9)=Delta_Bchl_ActSBPase;
% RDMB(10)=Delta_Bchl_ActPRK;
% RDMB(11)=Delta_Bchl_ActRubisco;

%vCO2s=reaction_flux(3);

% %leaf air space 20%
% Vol_airspace=0.04;%L
% Molar_Volume=22.4/273*(WeatherTemperature+273);
% LM_DYDT(m)=(vCO2s-Ameta)/(Vol_airspace/Molar_Volume);