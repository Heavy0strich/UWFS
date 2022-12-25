% MC Pressure Minimum Profile

% This script is used to compute the minimum Master Cylinder pressure
% profile needed to lock the wheels. This will be used to spec an actuator
% for the brake system. There are two computations done here. One is done
% for a static hydraulic system and the other utilizes the brake controller
% with high gains and high error. This is to ensure that the
% actuator/interface design is specced to ensure it is capable of going
% achieving the desired MC pressured with room for increasing if desired. 

% The computation: The greatest amount of torque to lock the wheels is seen
% at the highest vehicle speeds. This is due to aerodynamics and the help
% of downforce which helps with wheel grip. Therefore, the calculations
% will be done at the max speed and thus max torque. 


STAT = load('T32CarParamObj.mat','Brakes');
DYN1 = load('SimulationResults1Obj.mat','P_MC','k','lamb');
DYN2 = load('SimulationResults2Obj.mat','P_MC','k','lamb');
DYN3 = load('SimulationResults3Obj.mat','P_MC','k','lamb');
DYN4 = load('SimulationResults4Obj.mat','P_MC','k','lamb');
DYN5 = load('SimulationResults5Obj.mat','P_MC','k','lamb');

% Static system computation:
% in static system, Pw = Pmc
statPmc = (STAT.Brakes.LkTorque_F(1,end) / STAT.Brakes.Kb) + STAT.Brakes.Ppo;


% Dynamic system with high gain/error
dynPmc1 = max(DYN1.P_MC(1,:));
dynK1 = DYN1.k;
dynLamb1 = DYN1.lamb;


dynPmc2 = max(DYN2.P_MC(1,:));
dynK2 = DYN2.k;
dynLamb2 = DYN2.lamb;

dynPmc3 = max(DYN3.P_MC(1,:));
dynK3 = DYN3.k;
dynLamb3 = DYN3.lamb;


dynPmc4 = max(DYN4.P_MC(1,:));
dynK4 = DYN4.k;
dynLamb4 = DYN4.lamb;


dynPmc5 = max(DYN5.P_MC(1,:));
dynK5 = DYN5.k;
dynLamb5 = DYN5.lamb;



dynLambarray = [dynLamb1,dynLamb2,dynLamb3,dynLamb4,dynLamb5];
dynPmcarray = [dynPmc1,dynPmc2,dynPmc3,dynPmc4,dynPmc5];
%dynPmcarrayPSI = sort(dynPmcarray./6.895);
dynPmcarrayPSI = linspace(29,253.8,10);

figure()
yyaxis left
plot(dynLambarray,dynPmcarray)
set(gca,'YLim',[200, 1750])
xlabel('Lamda Gain')
ylabel('Pressure (kPa)')
title('Minimum MC Pressure to Lock Brakes vs. Lambda Gain')
hold on
plot(linspace(0,2.75,20),ones(1,20)*statPmc,'r-')
xlim([0.25 2.75])
legend('Dynamic','Static')

yyaxis right
set(gca,'YTick',dynPmcarrayPSI,'YLim',[29,253.8])
ytickformat('%.2f')
ylabel('Pressure (psi)')

yyaxis left
grid on



