%% Push Calorimeter Infusion Script - Pull Eq. Eqns.
% Katherine Ruud, MEI Research Ltd.
% This is the recommended way to analyze infusion data for the least amount
% of confusion. 
% This script is designed to calculate and simulate "Push" calorimeter eqns
% It calculates VO2 and VCO2 expected values based on the equivalent Pull
% calorimeter equations
% Outflow Rate = Inflow Rate + MFCs, Inflow Rate calc'd by Volume Haldane
% Outflow CO2 and O2 concentrations are calculated based on MFC values and
% a continous 'mixing tank' model for comparison to actual values.

% Input: Import CalRQ .csv data file as column vectors

% Output: Calculates gas concentration based on initial conditions and MFC 
% flows, plots comparisons to logged data

%% Constants

% InflowRate = InflowRate0;
CO2MFC = 2;
N2MFC = 3;

% Size of centered derivative (minutes)
dsize = 6;

% Corrections
Volume = 33100;

% Size of centered derivative (minutes)
dsize = 6;

% Correction Type
Type = "linear";
% Type = "gainoffset";

% chamber 1
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

% interpolation corrections
O2_in_interp = [
20.09495241
20.19401049
20.27787682
20.36540202
20.45474048
20.54640302
20.63871648
20.73496555
20.83503758
20.93136859
21.01106976
21.15095216
];
CO2_in_interp = [
0.000748668
0.096495943
0.19119915
0.289502001
0.390341035
0.492084718
0.596037341
0.700465933
0.804075924
0.906667964
1.008685726
];
O2_out_interp = [
20.09219475
20.18906286
20.27256644
20.3601434
20.44930248
20.54054355
20.63263749
20.72821319
20.82736094
20.9250207
21.00282749
21.14274324
];
CO2_out_interp = [
0.000650389
0.097367861
0.194177254
0.293466725
0.394144093
0.494886716
0.597404727
0.700300937
0.802221329
0.903496393
1.003704271

];
O2_ref = [
19.9998
20.1005
20.1994
20.2999
20.3997
20.5006
20.5994
20.6999
20.7997
20.9004
20.9992
21.1571
];
CO2_ref = [
0.0000
0.1000
0.2000
0.3000
0.4000
0.5000
0.6000
0.7000
0.8000
0.9000
1.0000
];

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
% [~,N2MFC] = max(mean(MFCMatrix));
% CO2Matrix = MFCMatrix;
% CO2Matrix(:,N2MFC) = 0;
% [~,CO2MFC] = max(mean(CO2Matrix));

% Find dt
dt = minutes(Time(2)-Time(1));

% Derivative Time
dtime = (dsize/2)*(1/dt);

% Apply Corrections and convert to fractional
if Type == "gainoffset"
    InflowO2Null = (InflowO2*O2Offset_in_m + O2Offset_in_b + O2Null)/100;
    OutflowO2Corrected = (OutflowO2*O2Offset_out_m + O2Offset_out_b)/100;
    InflowCO2Null = (InflowCO2*CO2Offset_in_m + CO2Offset_in_b + CO2Null)/100;
    OutflowCO2Corrected = (OutflowCO2*CO2Offset_out_m + CO2Offset_out_b)/100;
elseif Type == "linear"
    InflowO2Null = (interp1(O2_in_interp,O2_ref,InflowO2,'pchip') + O2Null)/100;
    InflowCO2Null = (interp1(CO2_in_interp,CO2_ref,InflowCO2,'pchip') + CO2Null)/100;
    OutflowO2Corrected = (interp1(O2_out_interp,O2_ref,OutflowO2,'pchip'))/100;
    OutflowCO2Corrected = (interp1(CO2_out_interp,CO2_ref,OutflowCO2,'pchip'))/100;
end

% % Filters
% windowSize = .5/dt;
% windowSize2 = 60/dt;
% windowSize3 = 3/dt;
% b = (1/windowSize)*ones(1,windowSize);
% c = (1/windowSize2)*ones(1,windowSize2);
% d = (1/windowSize3)*ones(1,windowSize3);
% a = 1;
% InflowO2Null = filter(c,a,InflowO2Null);
% OutflowO2Corrected = filter(d,a,OutflowO2Corrected);
% InflowCO2Null = filter(b,a,InflowCO2Null);
% OutflowCO2Corrected = filter(b,a,OutflowCO2Corrected);

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
    % (expected change based on MFCs, concentration and flows)
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

    % Filter O2 D
    % OutflowO2D = filter(b,a,OutflowO2D);
    
% Calculate Expected Values
VO2Exp = (InflowO2Null .* MFCMatrix(:,N2MFC)) ./ (InflowN2Null);
VCO2Exp = (MFCMatrix(:,CO2MFC) .* (InflowN2Null) - InflowCO2Null .* MFCMatrix(:,N2MFC)) ./ (InflowN2Null);

% Calculate "Haldane" Outflow Rate
HInflow = (OutflowRate .* OutflowN2Corrected + OutflowN2D * Volume) ./ InflowN2Null;

% Calculate VO2 VCO2 RQ
VO2 = HInflow .* InflowO2Null - OutflowRate .* OutflowO2Corrected - OutflowO2D * Volume;
VCO2 = -( (HInflow .* InflowCO2Null - OutflowRate .* OutflowCO2Corrected) - OutflowCO2D * Volume);
RQ = VCO2 ./ VO2;

% VCO2test = -( (HInflow .* InflowCO2Null - OutflowRate .* OutflowCO2) - OutflowCO2DSim * Volume);

% Calculate VO2 VCO2 Error
VO2E = abs(((VO2 - VO2Exp) ./ (VO2Exp)) * 100);
VCO2E = abs(((VCO2 - VCO2Exp) ./ (VCO2Exp)) * 100);

%% Plot Data
close all
% plot(medfilt1(HInflow,500) .* medfilt1(InflowO2Null,500) - OutflowRate .* medfilt1(OutflowO2Corrected,500) - medfilt1(OutflowO2D,500) * Volume)
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
plot(Time,OutflowO2D,Time,OutflowO2DSim,Time,OutflowO2DSim)
legend('Outflow D','Simulated Outflow D','Simulated D')
title('O2')

figure
plot(Time,OutflowCO2D,Time,OutflowCO2DSim)
legend('Outflow D','Simulated Outflow D','Simulated D')
title('CO2')

figure
plot(Time,OutflowN2D,Time,OutflowN2DSim,Time,OutflowN2DSim)
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
 
% Plot VO2/VCO2 Error
figure
plot(Time,VO2E,Time,ones(length(Time),1)*4)
legend('VO2 Error','4% Error')
ylim([0,100]);
title('VO2 Error')

figure
plot(Time,VCO2E,Time,ones(length(Time),1)*4)
legend('VCO2 Error','4% Error')
ylim([0,100]);
title('VCO2 Error')


