% Retime a CalRQ data file to minute-by-minute data

[GFileName,PathName] = uigetfile('*csv','File to convert','MultiSelect','on');

% Handle single/multi file select
if PathName == 0
    disp('Error: No file selected')
    return
elseif ischar(GFileName)
    VarNumMovs = 1;
    MFileName{1,1} = GFileName;
else
    VarNumMovs = length(GFileName);
    MFileName = GFileName;
end

for i=1:VarNumMovs
    %% Analyze code
    FileName=MFileName{i};
    
%% Initialize variables.
filename = [PathName,FileName];
delimiter = ',';
startRow = 4;

%% Format for each line of text:
%   column1: datetimes (%{MM/dd/yy HH:mm:ss}D)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: text (%q)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: double (%f)
%	column24: double (%f)
%   column25: double (%f)
%	column26: double (%f)
%   column27: double (%f)
%	column28: double (%f)
%   column29: double (%f)
%	column30: double (%f)
%   column31: text (%q)
%	column32: text (%q)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%{MM/dd/yy HH:mm:ss}D%f%f%f%f%q%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
T = table(dataArray{1:end-1}, 'VariableNames', {'Time','VO2','VCO2','RQ','MR','Activity','InflowO2','InflowCO2','InflowH2O','OutflowO2','OutflowCO2','OutflowH2O','InflowRate','OutflowRate','ChamberTemperature','ChamberPressure','ChamberHumidity','InflowH2O1','LabTemp','LabHumidity','BIOSFlowL','BIOSTempL','BIOSBPL','BIOSFlowH','BIOSTempH','BIOSBPH','MFCFlow_1','MFCFlow_2','MFCFlow_3','MFCFlow_4','Alarms','Notes'});

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% SCDelayTestLow45LPM1Infusion170911.Time=datenum(SCDelayTestLow45LPM1Infusion170911.Time);

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

% T = table(Time,InflowO2,InflowCO2,OutflowO2,OutflowCO2,InflowRate,ChamberPressure,ChamberTemperature,ChamberHumidity,OutflowH2O,InflowH2O,MFCFlow_1,MFCFlow_2,MFCFlow_3,MFCFlow_4,BIOSFlowL,BIOSFlowH);
% T.Properties.VariableNames = {'Time' 'O2In' 'CO2In' 'O2Out' 'CO2Out' 'InflowRate' 'ChamberPressure' 'ChamberTemp' 'ChamberHumidity'  'OutflowH2O' 'InflowH2O' 'MFC1' 'MFC2' 'MFC3' 'MFC4' 'BIOSFlowL' 'BIOSFlowH'};
T = table2timetable(DelayTestLow11Infusion171026);
T = retime(T,'minutely');
T = timetable2table(T);
writetable(T,'DelayTestLow11Infusion171026.xlsx');

end