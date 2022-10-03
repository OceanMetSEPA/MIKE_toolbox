function [ op ] = logInfo( logFile )
% Provide information about MIKE model run (extracted from its log file)
%
% INPUT: either
% a) filename of log file
% b) directory containing model run
%
% OUTPUT:
% struct containing information about when model was run, how long it took,
% and whether it completed properly
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   logInfo.m  $
% $Revision:   1.1  $
% $Author:   ted.schlicke  $
% $Date:   Feb 02 2018 11:36:18  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin~=1
    help MIKE.logInfo
    return
end

if isdir(logFile)
    logFile=fileFinder(logFile,'.log');
end
if iscell(logFile)
    logFile=char(logFile);
end
if min(size(logFile))>1
    error('Multiple log files found! Please specify unique file')
end
txt=readTxtFile(logFile);
logDate=strtrim(txt{1});
logTime=strtrim(txt{2});
logDateTime=datenum(strcat(logDate,logTime),'yyyymmddHHMMSS');

%
startDateRow=stringFinder(txt,'Total','output','index');
str=strsplit(txt{startDateRow(end)},' ');
runTime=str2double(str{end})/3600; % number of hours
%degreeConverter(runTime,3)

%fprintf('Model run time = %s\n',days2String(runTime/24))

modelPath=logFile(1:max(regexp(logFile,'\')));
modelPath=driveLetter2Hostname(modelPath);
fprintf('Model Path = ''%s''\n',modelPath)
di=dir(logFile);

% Memory type:
isOpenMP=stringFinder(txt,'Computer','output','any');
if isOpenMP
    memoryType='OpenMP';
else
    memoryType='MPI';
end

op=struct;
op.Path=char(modelPath);
op.StartTime=datestr(logDateTime);
op.EndTime=datestr(di.datenum);
op.runTime=days2String(runTime/24);
op.memoryType=memoryType;
op.txt=txt;

if stringFinder(txt{end},{'Abnormal','cancelled'},'type','or','output','bool')
    warning('OH DEAR; run finished abnormally')
    op.Success=false;
else
    op.Success=true;
end
