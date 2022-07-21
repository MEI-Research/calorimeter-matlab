%% Push Calorimeter Infusion Script - Corrected Push Method
% Katherine Ruud, MEI Research Ltd.
% This script is designed to calculate and simulate "Push" calorimeter eqns
% It calculates VO2 and VCO2 expected values based based on Push
% calorimeter equations.
% VO2 and VCO2 derivative terms are corrected so that they remain constant
% for constant blender MFC values.
% Outflow CO2 and O2 concentrations are calculated based on MFC values and
% a continous 'mixing tank' model

% Input: Import CalRQ .csv data file as column vectors

% Output: Calculates gas concentration based on initial conditions and MFC 
% flows, plots comparisons to logged data

%% Constants
% Corrections
Volume = 30000;
O2Offset_in_m = 1;
O2Offset_in_b = 0;
O2Offset_out_m = 1;
O2Offset_out_b = 0;
O2Null = 0;

CO2Offset_in_m = 1;
CO2Offset_in_b = 0;
CO2Offset_out_m = 1;
CO2Offset_out_b = 0;
CO2Null = 0;

% Size of centered derivative (minutes)
dsize = 8;

% MFC assignments automatically set based on avg flow

%% Calculations
% Clear old variables that may cause issues
clear OutflowO2D OutflowCO2D OutflowN2D

% Replace NaN
% Could also interpolate missing values but at 1 minute intervals we
% probably won't see NaN from any instruments
MFCFlow_1(isnan(MFCFlow_1)) = 0;
MFCFlow_2(isnan(MFCFlow_2)) = 0;
MFCFlow_3(isnan(MFCFlow_3)) = 0;
MFCFlow_4(isnan(MFCFlow_4)) = 0;
InflowRate(isnan(InflowRate)) = 0;

% Assume N2 MFC is Max
% Can also manually assign MFCs if you are doing a profile with only CO2 or
% extreme RQs (VCO2 >> VO2)
MFCMatrix = [MFCFlow_1,MFCFlow_2,MFCFlow_3,MFCFlow_4];
[~,N2MFC] = max(mean(MFCMatrix));
CO2Matrix = MFCMatrix;
CO2Matrix(:,N2MFC) = 0;
[~,CO2MFC] = max(mean(CO2Matrix));

% Find dt
dt = minutes(Time(2)-Time(1));

% Derivative Time
dtime = (dsize/2)*(1/dt);

% Apply Corrections and convert to fractional
InflowO2Null = (InflowO2*O2Offset_in_m + O2Offset_in_b + O2Null)/100;
OutflowO2Corrected = (OutflowO2*O2Offset_out_m + O2Offset_out_b)/100;
InflowCO2Null = (InflowCO2*CO2Offset_in_m + CO2Offset_in_b + CO2Null)/100;
OutflowCO2Corrected = (OutflowCO2*CO2Offset_out_m + CO2Offset_out_b)/100;

% Initialize Variables
OutflowO2Sim = zeros(length(InflowO2)+1,1);
OutflowCO2Sim = zeros(length(InflowO2)+1,1);

% Calculate Outflow Rate
OutflowRate = InflowRate + MFCMatrix(:,N2MFC) + MFCMatrix(:,CO2MFC);
OutflowN2Corrected = 1 - OutflowO2Corrected - OutflowCO2Corrected;
InflowN2Null = 1 - InflowO2Null - InflowCO2Null;

% Model Chamber
for i = 1:length(InflowO2Null)
    % What are our initial conditions?
    if i == 1
        OutflowO2Sim(1) = OutflowO2Corrected(1);
        OutflowCO2Sim(1) = OutflowCO2Corrected(1);
        OutflowN2Sim(1) = OutflowN2Corrected(1);
    end
    % Fix for divide by zero
    if OutflowRate(i) == 0
        OutflowRate(i) = 1;
    end
    
    % Calculate Model Constant
    ConO2 = - OutflowO2Sim(i) * OutflowRate(i) + InflowRate(i) * InflowO2Null(i);
    ConCO2 = - OutflowCO2Sim(i) * OutflowRate(i) + InflowRate(i) * InflowCO2Null(i) + MFCMatrix(i,CO2MFC);
    ConN2 = - OutflowN2Sim(i) * OutflowRate(i) + InflowRate(i) * InflowN2Null(i) + MFCMatrix(i,N2MFC);
    
    % Calculate Next Value
    OutflowO2Sim(i+1,1) = ( 1 / OutflowRate(i) ) * (- ConO2 * exp( -OutflowRate(i) * dt / Volume ) + InflowRate(i) * InflowO2Null(i));
    OutflowCO2Sim(i+1,1) = ( 1 / OutflowRate(i) ) * (- ConCO2 * exp( -OutflowRate(i) * dt / Volume ) + InflowRate(i) * InflowCO2Null(i) + MFCMatrix(i,CO2MFC));
    OutflowN2Sim(i+1,1) = ( 1 / OutflowRate(i) ) * (- ConN2 * exp( -OutflowRate(i) * dt / Volume ) + InflowRate(i) * InflowN2Null(i) + MFCMatrix(i,N2MFC));
      
    % Calculate Simulation Derivatives
    OutflowO2DSim(i,1) = (InflowRate(i) * InflowO2Null(i) - OutflowO2Sim(i) * OutflowRate(i)) / Volume;  
    OutflowCO2DSim(i,1) = (InflowRate(i) * InflowCO2Null(i) + MFCMatrix(i,CO2MFC) - OutflowCO2Sim(i) * OutflowRate(i)) / Volume;    
    OutflowN2DSim(i,1) = (InflowRate(i) * InflowN2Null(i) + MFCMatrix(i,N2MFC) - OutflowN2Sim(i) * OutflowRate(i)) / Volume;
   
    % Calculate Real Derivative
    if i > dtime && i < length(OutflowO2Corrected)-dtime
        X = [ones(2*dtime+1,1) transpose(1:(2*dtime+1))*dt];
        
        OutO2D = X\OutflowO2Corrected(i-dtime:i+dtime);
        OutCO2D = X\OutflowCO2Corrected(i-dtime:i+dtime);
        OutN2D = X\OutflowN2Corrected(i-dtime:i+dtime);
        
        OutflowO2D(i,1) = OutO2D(2);
        OutflowCO2D(i,1) = OutCO2D(2);
        OutflowN2D(i,1) = OutN2D(2);
    else
        OutflowO2D(i,1) = 0;
        OutflowCO2D(i,1) = 0;
        OutflowN2D(i,1) = 0;
    end
end

    OutflowO2Sim = OutflowO2Sim(1:end-1);
    OutflowCO2Sim = OutflowCO2Sim(1:end-1);
    OutflowN2Sim = OutflowN2Sim(1:end-1);

% Calculate Expected Values
VO2Exp = (InflowRate .* InflowO2Null .* MFCMatrix(:,N2MFC)) ./ ( MFCMatrix(:,N2MFC) + InflowRate .* InflowN2Null);
VCO2Exp = InflowRate .* ( InflowCO2Null .* MFCMatrix(:,N2MFC) - MFCMatrix(:,CO2MFC) .* InflowN2Null) ./ ( - MFCMatrix(:,N2MFC) - InflowRate .* InflowN2Null);

% Calculate Derivative Correction
DCorr = (InflowRate - VO2Exp + VCO2Exp) ./ (InflowRate + MFCMatrix(:,N2MFC) + MFCMatrix(:,CO2MFC));

% Calculate "Haldane" Outflow Rate
HOutflow = (InflowRate .* InflowN2Null - OutflowN2D .* DCorr * Volume) ./ OutflowN2Corrected;

% Calculate VO2 VCO2 RQ
VO2 = InflowRate .* InflowO2Null - HOutflow .* OutflowO2Corrected - OutflowO2D .* DCorr * Volume;
VCO2 = -( (InflowRate .* InflowCO2Null - HOutflow .* OutflowCO2Corrected) - OutflowCO2D .* DCorr * Volume);
RQ = VCO2./VO2;

% Calculate VO2 VCO2 Error
VO2E = abs(((VO2 - VO2Exp) ./ (VO2Exp)) * 100);
VCO2E = abs(((VCO2 - VCO2Exp) ./ (VCO2Exp)) * 100);

%% Plot Data
close all

% Plot Simulation
figure
plot(Time,OutflowO2Sim,Time,OutflowO2Corrected)
legend('Simulated','Real')
title('O2')

figure
plot(Time,OutflowCO2Sim,Time,OutflowCO2Corrected)
legend('Simulated','Real')
title('CO2')

%%figure
%plot(Time,OutflowN2Sim,Time,OutflowN2Corrected)
%legend('Simulated','Real')
%title('N2')

% Plot Derivatives
figure
plot(Time,OutflowO2D,Time,OutflowO2DSim)
legend('Outflow D','Simulated Outflow D','Simulated D')
title('O2')

figure
plot(Time,OutflowCO2D)
% plot(Time,OutflowCO2D,Time,OutflowCO2DSim,Time,OutflowCO2DSim2)
legend('Outflow D','Simulated Outflow D','Simulated D')
title('CO2')

figure
plot(Time,OutflowN2D,Time,OutflowN2DSim)
legend('Outflow D','Simulated Outflow D','Simulated D')
title('N2')

% Plot VO2/VCO2 Data
figure
plot(Time,VO2,Time,VO2Exp)
legend('VO2','VO2 Exp')
title('VO2')

figure
plot(Time,VCO2,Time,VCO2Exp)
legend('VCO2','VCO2 Exp')
title('VCO2')
% 
% % Plot VO2/VCO2 Error
figure
plot(Time,VO2E,Time,ones(length(Time),1)*4)
legend('VO2','VO2 Exp')
ylim([0,100]);
title('VO2')

figure
plot(Time,VCO2E,Time,ones(length(Time),1)*4)
legend('VCO2','VCO2 Exp')
ylim([0,100]);
title('VCO2')

