%Sensitivity analysis of photosynthetic parameters
clear all;
SpeciesCol=1; %1Maize	2sorghum	3sugarcane
for i=1:19
    Rp1=RAC4leafMetaDrive(i,1,SpeciesCol);
    Rp2=RAC4leafMetaDrive(i,-1,SpeciesCol);
    Rp0=RAC4leafMetaDrive(i,0,SpeciesCol);
    CC=(Rp1(:,2)-Rp2(:,2))./Rp0(:,2)/0.02;
    expr=['CC' num2str(i) '= CC;'];
    eval(expr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output: Sensitivity coeffcient of photosynthetic parameters
% ki	1/TaoRubisco	[RP]	Rubisco	Vcmax	MDH	PEPC	1/Taopepc	PPDK	MDH	ME	DAPDH	SBPase	FBPase	PRK	I	Mutase&Enolase	Jmax Gs A(CO2 uptake rate)		
RCC=[CC1,CC2,CC3,CC5,CC4,CC6,CC7,CC8,CC9,CC10,CC11,CC12,CC13,CC14,CC15,CC16,CC17,CC18,CC19,Rp0(:,2)]; 
Timex=Rp1(:,1);
%Generate figures
SensitivityFig(Timex,RCC);
SensitivityRCC=[Timex,RCC];% Add Column 1 (time)
dlmwrite('../Results/SensitivityRCC.txt', SensitivityRCC, '\t')