% Format result table into Excel sheet

T = table(Time,VO2*1000,VCO2*1000,VO2Exp*1000,VCO2Exp*1000,VO2E,VCO2E,InflowO2Null*100,InflowCO2Null*100,OutflowO2Corrected*100,OutflowCO2Corrected*100,InflowRate, ChamberPressure, MFCFlow_1, MFCFlow_2, MFCFlow_3, MFCFlow_4);
T.Properties.VariableNames = {'Time' 'VO2' 'VCO2' 'VO2Expected' 'VCO2Expected' 'VO2Error' 'VCO2Error' 'InflowO2' 'InflowCO2' 'OutflowO2' 'OutflowCO2' 'InflowRate' 'ChamberPressure' 'MFC1' 'MFC2' 'MFC3' 'MFC4'};
T = table2timetable(T);
T = retime(T,'minutely');
T = timetable2table(T);
writetable(T,'data4.xls');