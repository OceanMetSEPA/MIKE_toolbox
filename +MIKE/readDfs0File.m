function [dfs0Output ] = readDfs0File( dfs0FileName,varargin )
% Extract data from MIKE Dfs0 file
%
% This function extracts time-series data from a MIKE Dfs0 file
% (corresponding to a single point)
%
% INPUT:
% dfs0FileName - file used by MIKE with extension .dfs0
%
% Optional Inputs:
% info (false) - return information about file rather than data (e.g.
%                data parameters, time period, location names)
% parameters [] - return specified parameters
% locations []  - return data for specified location
% combineDataAtSingleLocation [true] - group data for single parameter together
% dataAsFields [false] - store data as struct field rather than separate struct
%
% OUTPUT:
% struct containing metadata about Dfs0 file and data in various formats,
% depending on optional inputs above. Let the number of locations be Nloc
% and the number of parameters by Nparam.
%
% Case A:       dataAsFields == false %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data stored as struct array within 'Data' field:
%
% 1) if combineDataAtSingleLocation == false:
%       'Data' is a struct array of Nparam * Nloc structs. Each struct
%       correponds to a separate location and parameter. In each struct,
%       the data is contained within the 'Value' field.
%
% 2) if combineDataAtSingleLocation == true:
%       'Data' is a struct array of Nloc structs. Each struct corresponds
%       to the data at a given location. The location is specified in the
%       'Location' field, and data for individual parameters are stored in
%       fields with format 'Parameter(unit)'.
%
% Case B:       dataAsFields == true %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) if combineDataAtSingleLocation == false:
%       Data stored in Nlocs * Nparam separate fields, with fieldname of
%       format 'Location_Parameter(Unit)'
%
% 2) if combineDataAtSingleLocation == true:
%       Data stored in Nloc fields, each with fieldname given by location.
%       Each of these fields consists of a struct where the data is stored
%       in field with fieldname of format 'Parameter(unit)'.
%
% Note 1 - there may be too many output options here but i wasn't sure which
% would be most useful! Might be able to simplify in future.
%
% Note 2 - Matlab struct fieldnames can't use '(' or ')', so fieldnames
% have been passed through genvarname function. This can be reversed using
% PVCS function 'ungenvarname'. It was deemed useful to retain the units
% for plotting purposes etc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   readDfs0File.m  $
% $Revision:   1.2  $
% $Author:   ted.schlicke  $
% $Date:   Feb 02 2018 11:40:10  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0
    help MIKE.readDfs0File
    return
end

options=struct;
options.info=false;
options.parameters=[];
options.locations=[];
options.combineDataAtSingleLocation=1;
options.dataAsFields=0;
options=checkArguments(options,varargin);

if ~exist(dfs0FileName,'file')
    error('File ''%s'' not found',dfs0FileName)
end

if ~stringFinder(dfs0FileName,'.dfs0','type','end','output','bool')
    error('File ''%s'' not a .dfs0 file',dfs0FileName)
end

timeScaleFactor=24*60*60;

% A script loading a dfs0 file and plotting the time series, using the
% Dfs0Util for reading the dfs0 data efficiently. Without the use of the
% Dfs0Util, performance is very bad.

%dfs0FileName='C:\MIKEZeroProjectsSEPA\ShunaSound\2003\PART_IncludeHDOutput\sliceDischarge.dfs0';

if ~stringFinder(dfs0FileName,'.dfs0','type','end','output','bool')
    error('Not a dfs0 file! Aborting...')
end
% 
NETaddAssembly('DHI.Generic.MikeZero.DFS.dll');
import DHI.Generic.MikeZero.DFS.*;
import DHI.Generic.MikeZero.DFS.dfs0.*;
% OLD STUFF
% NET.addAssembly('DHI.Generic.MikeZero.DFS');
% import DHI.Generic.MikeZero.DFS.*;
% import DHI.Generic.MikeZero.DFS.dfs0.*;
% FROM readDfsuFile- doesn't work! 
% NETaddAssembly('DHI.Generic.MikeZero.DFS.dll');
% NETaddAssembly('DHI.Generic.MikeZero.DFS.dll');
% NETaddAssembly('DHI.Generic.MikeZero.EUM.dll');
% import DHI.Generic.MikeZero.DFS.*;
% import DHI.Generic.MikeZero.*.*;


dfs0File  = DfsFileFactory.DfsGenericOpen(dfs0FileName);

%% Dfs0 Items
NItems=dfs0File.ItemInfo.Count;
dfs0Items=arrayfun(@(i)char(dfs0File.ItemInfo.Item(i-1).Name),1:NItems,'Unif',0)';
% Data with location info seem to have format 'Location: Parameter'
% Data without location info have format 'Parameter'
% Split string into these two things:
% Locations & Parameters
str=regexp(dfs0Items,':','split');
str=vertcat(str{:});
if size(str,2)==1 % No location info
    dfs0Locations=[];
    dfs0Parameters=strtrim(unique(str));
else
    dfs0Locations=unique(str(:,1));
    dfs0Parameters=strtrim(unique(str(:,2)));
end

%% Read times and data for all items
% Use the Dfs0Util for bulk-reading all data and timesteps
dd = double(Dfs0Util.ReadDfs0DataDouble(dfs0File));

% Time stuff -
timeColumn = dd(:,1); % first column of data - seem to be in seconds past start time
% In dfs0 file generated from Mike -> File -> New, timestep seems to be
% minutes, even though dfs0File.FileInfo.TimeAxis.TimeUnit = 'eumUsec' ?!
t0=dfs0File.FileInfo.TimeAxis.StartDateTime;
y=t0.Year;
m=t0.Month;
d=t0.Day;
H=t0.Hour;
M=t0.Minute;
S=t0.Second;
vec=double([y,m,d,H,M,S]);
t0=datenum(vec);
% Use above parameters to generate timings of data (in datenums)
dateTime=t0+timeColumn/timeScaleFactor;

%% Prepare struct with basic dfs0 information
metaData=struct('FileName',dfs0FileName,'FileTitle',char(dfs0File.FileInfo.FileTitle));
metaData.Locations=dfs0Locations;
metaData.Parameters=dfs0Parameters;
metaData.Items=dfs0Items;

dfs0Output=struct('MetaData',metaData);

% Does user just want info about file?
if options.info
    fprintf('Returning info only\n')
    return
end

%% Check data requirements
if ~isempty(options.parameters)
    fprintf('Checking parameters...\n')
    kParam=stringFinder(dfs0Items,options.parameters,'output','index','type','or','ignorecase',1);
    if isempty(kParam)
        fprintf('Valid parameters:\n')
        disp(dfs0Parameters)
        error('Specified parameter not found; please select one of the above')
    end
else
    kParam=1:NItems; % All items
end

if ~isempty(options.locations)
    fprintf('Checking locations...\n')
    kLoc=stringFinder(dfs0Items,options.locations,'output','index','type','or','ignorecase',1);
    if isempty(kLoc)
        if ~isempty(dfs0Locations)
            fprintf('Valid locations:\n')
            disp(dfs0Locations)
        else
            fprintf('No location information available; please don''t filter by location\n')
        end
        error('Specified location not found')
    end
else
    kLoc=1:NItems;
end

itemIndices2Extract=intersect(kParam,kLoc);

pointData=dd(:,2:end);

%% Load each location and parameter's data into struct array
NItems=length(itemIndices2Extract);
dataStructArray=cell(NItems,1); % Convert to struct array after loop
for itemIndex=1:NItems
    ithItemIndex=itemIndices2Extract(itemIndex);
    %    fprintf('Loading item %d of %d...\n',i,NItems)
    item=dfs0File.ItemInfo.Item(ithItemIndex-1);
    name=char(item.Name);  % name given by MIKE, for example "Shuna_Castle_Bay_2000: Surface elevation"
    %        unit=char(item.Quantity.Unit);
    splitName=strsplit(name,':'); % First bit of name (Observation point name)
    if length(splitName)==1
        locationName=[];
        paramName=strtrim(splitName{1});
    else
        locationName=strtrim(splitName{1});
        paramName=strtrim(splitName{end});
    end
    quantity=item.Quantity;
    unitDescription=char(quantity.UnitDescription);
    unitAbbreviation=char(quantity.UnitAbbreviation);
    %    fieldName=genvarname(sprintf('%s(%s)',paramName,unitAbbreviation));
    s=struct('Dfs0Name',name,'Location',locationName,'Parameter',paramName,'DateTime',dateTime,'Value',pointData(:,ithItemIndex),'UnitDescription',unitDescription,'UnitAbbreviation',unitAbbreviation);
    dataStructArray{ithItemIndex}=s;
end
dataStructArray=vertcat(dataStructArray{:}); %cell -> struct array

%% bundle structs in corresponding to same observation point together
if options.combineDataAtSingleLocation
    if isempty([dataStructArray.Location])
        %        warning('No location information available; can''t combine data')
        data=dataStructArray;
    else
        %   fprintf('Combine data by location\n')
        locs={dataStructArray.Location};
        locationNames=unique(locs);
        Nop=length(locationNames);
        dataStructByLocation=cell(length(locationNames),1);
        for locationIndex=1:Nop
            locationName=locationNames{locationIndex};
%            fprintf('%d of %d : %s\n',locationIndex,Nop,locationName)
            si=dataStructArray(strcmp(locationName,locs));
            s=struct('Location',locationName,'DateTime',dateTime);
            for parameterIndex=1:length(si)
                sp=si(parameterIndex);
                param=sprintf('%s(%s)',sp.Parameter,sp.UnitAbbreviation);
                fn=genvarname(param);
                s.(fn)=sp.Value;
            end
            dataStructByLocation{locationIndex}=s;
        end
        try % to bundle cell array to struct array
            dataStructByLocation=vertcat(dataStructByLocation{:});
            data=dataStructByLocation;
        catch
            warning('Failed to bundle structs together')
            % check fieldnames
            Nd=length(dataStructByLocation);
            fn=cellfun(@(x)ungenvarname(fieldnames(x)),dataStructByLocation,'Unif',0);
            fnTally=tally(vertcat(fn{:}));
            k=[fnTally{2,:}]~=Nd;
            fnTally=fnTally(:,k);
            if ~isempty(fnTally)
                warning('Mismatched fieldnames!')
                disp(fnTally)
            end
            fprintf('Try to use fieldnames from 1st struct for all others...\n')
            fn1=fieldnames(dataStructByLocation{1});
            try
                for dataIndex=2:Nd
                    dataStructByLocation{dataIndex}=renameStructFields(dataStructByLocation{dataIndex},fn1);
                end                
                dataStructByLocation=vertcat(dataStructByLocation{:});
                data=dataStructByLocation;
                fprintf('Success!\n')
            catch
                warning('Still can''t create struct array...')
                data=dataStructByLocation;
            end
        end
    end
else
    data=dataStructArray;
end

%%
if options.dataAsFields
    Ns=length(data);
    if options.combineDataAtSingleLocation
        struct2Add=struct;
    else
        struct2Add=struct('DateTime',dateTime);
    end
    for structIndex=1:Ns
        ithDataStruct=data(structIndex);
        if options.combineDataAtSingleLocation
            fni=MIKE.getFieldName(ithDataStruct.Location);
            dfs0Output.(fni)=rmfield(ithDataStruct,'Location');
        else
            fn=genvarname(sprintf('%s_%s(%s)',getMikeFieldName(ithDataStruct.Location),ithDataStruct.Parameter,ithDataStruct.UnitAbbreviation));
            struct2Add.(fn)=ithDataStruct.Value;
        end
    end
    dfs0Output=mergestruct(dfs0Output,struct2Add);
else
    dfs0Output.Data=data;
end

dfs0File.Close()

