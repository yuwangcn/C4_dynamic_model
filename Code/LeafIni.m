function con = LeafIni()
global WeatherTemperature;
global WeatherRH;
global Air_CO2;
global BallBerryInterceptC4;
Tleaf=WeatherTemperature;%energy balance 
H2Oou=0;
CO2in=0;
Ci = 0.7 * Air_CO2;%Initial Ci u moles/mole solublility check
Gs = BallBerryInterceptC4; % Initial stomatal conductance moles/m2 leaf area/s
Cb= 0.9* Air_CO2;
ESaturation = 0.611 * exp(17.502 * WeatherTemperature / (WeatherTemperature + 240.97));%KPa 
Ea = WeatherRH * ESaturation;%KPa
Eb=Ea;
con=zeros(1,7);
con(1)=Ci;
con(2)=Cb;
con(3)=Eb;
con(4)=Gs;
con(5)=Tleaf;
con(6)=H2Oou;
con(7)=CO2in;
