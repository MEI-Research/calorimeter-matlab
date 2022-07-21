%% Push Calorimeter Human Script
% Katherine Ruud, MEI Research Ltd.
% This script calculates VO2 VCO2 RQ for human participants in a 'Push'
% calorimeter

% Input: Import CalRQ .csv data file as column vectors

% Output: Calculates VO2 VCO2 RQ based on equivelent equations to 
% Brown 1984 paper
tic 

%% Constants
% Corrections
Volume = 5000;

% Size of centered derivative (minutes)
dsize = 3;

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
20.05645041
20.15636869
20.25700099
20.3588109
20.45792417
20.55738266
20.66243646
20.76184585
20.86770205
20.96425726
21.06614196
21.2
];
CO2_in_interp = [
0.001929841
0.097271417
0.191591996
0.289089659
0.389456217
0.490423605
0.594143218
0.698076097
0.801219298
0.903467864
1.005660399
];
O2_out_interp = [
20.01633516
20.11763883
20.22048447
20.32282653
20.42408886
20.52467684
20.63055866
20.73096512
20.83811693
20.93625239
21.03962068
21.2
];
CO2_out_interp = [
-0.000630894
0.098487434
0.197554627
0.298733603
0.401476052
0.503985191
0.608734643
0.713453892
0.817358452
0.920496252
1.023233064
];
O2_ref = [
20.0000495
20.10002881
20.20060317
20.29949121
20.40026744
20.50036742
20.59912133
20.69990992
20.80051036
20.90034381
20.9993296
21.13218738
];
CO2_ref = [
0
0.100011643
0.199983796
0.300003805
0.400004039
0.49998422
0.600000061
0.699958156
0.799981122
0.89998216
0.999921918
];

%% Fix Excel time merge
if Time{1,1} == Time{2,1}
    TimeOld = Time;
    TimeDiff = (Time(60)-Time(1))/60;
    for i = 1:length(Time)
        TimeNew(i) = TimeOld(1) + TimeDiff*i - TimeDiff;
    end
    Time = transpose(TimeNew);
end        

%% Calculations
% Clear old variables that may cause issues
clear OutflowO2D OutflowCO2D OutflowN2D

% Replace NaN
% Could also interpolate missing values but at 1 minute intervals we
% probably won't see NaN from any instruments
InflowRate(isnan(InflowRate)) = 0;

% Find dt
dt = minutes(Time{2,1}-Time{1,1});

% Derivative Time
dtime = (dsize/2)*(1/dt);

% Apply Corrections and convert to fractional
if Type == "gainoffset"
    InflowO2Null = (InflowO2*O2Offset_in_m + O2Offset_in_b + O2Null)/100;
    OutflowO2Corrected = (OutflowO2*O2Offset_out_m + O2Offset_out_b)/100;
    InflowCO2Null = (InflowCO2*CO2Offset_in_m + CO2Offset_in_b + CO2Null)/100;
    OutflowCO2Corrected = (OutflowCO2*CO2Offset_out_m + CO2Offset_out_b)/100;
elseif Type == "linear"
    InflowO2Null = (interp1(O2_in_interp,O2_ref,InflowO2,'makima') + O2Null)/100;
    InflowCO2Null = (interp1(CO2_in_interp,CO2_ref,InflowCO2,'makima') + CO2Null)/100;
    OutflowO2Corrected = (interp1(O2_out_interp,O2_ref,OutflowO2,'makima'))/100;
    OutflowCO2Corrected = (interp1(CO2_out_interp,CO2_ref,OutflowCO2,'makima'))/100;
end
OutflowN2Corrected = 1 - OutflowO2Corrected - OutflowCO2Corrected;
InflowN2Null = 1 - InflowO2Null - InflowCO2Null;

% Derivative
for i = 1:length(InflowO2Null)
    
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

% Calculate "Haldane" Outflow Rate
HOutflow = (InflowRate .* InflowN2Null - OutflowN2D * Volume) ./ OutflowN2Corrected;

% Calculate VO2 VCO2 RQ
VO2 = InflowRate .* InflowO2Null - HOutflow .* OutflowO2Corrected - OutflowO2D * Volume;
VCO2 = -( (InflowRate .* InflowCO2Null - HOutflow .* OutflowCO2Corrected) - OutflowCO2D * Volume);
RQ = VCO2./VO2;

%% Plot Data
close all

figure
plot(Time{1:end,1},VO2)
legend('VO2','VO2 Exp')
% ylim([0,100]);
title('VO2')

figure
plot(Time{1:end,1},VCO2)
legend('VCO2')
% ylim([0,100]);
title('VCO2')

results = table(Time{1:end,1},InflowO2Null,InflowCO2Null,OutflowO2Corrected,OutflowCO2Corrected,InflowRate,VO2,VCO2);
%% Convert time interval
% test = table(Time,InflowO2Null,InflowCO2Null,OutflowO2Corrected,OutflowCO2Corrected,InflowRate)
% table2timetable(test) 
% test2 = retime(test,'minutely');

toc