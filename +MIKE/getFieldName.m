function [ outputStrings ] = getFieldName( inputStrings )
% Generate valid fieldname from MIKE item names
%
% INPUT:
% string (char or cellstr)
%
% OUTPUT:
% outputStrings - string (char) or strings (cellstr) suitable for struct fieldnames
%
% Function used by readDfs0File, readDfsuFile
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   getFieldName.m  $
% $Revision:   1.0  $
% $Author:   ted.schlicke  $
% $Date:   Jul 29 2016 11:22:02  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0
    help MIKE.getFieldName
    return
end

% Make sure we're working with cellstr:
if ischar(inputStrings)
    inputStrings=cellstr(inputStrings);
end

NStrings=length(inputStrings);
outputStrings=cell(NStrings,1);
for stringIndex=1:NStrings
    inputString=inputStrings{stringIndex};
    Nc=length(inputString);
    % Keep first letter, and capitalise letters following characters which
    % aren't letters
    for i=2:Nc
        if isletter(inputString(i)) && ~isletter(inputString(i-1))
            inputString(i)=upper(inputString(i));
        end
    end
    fn=uint8(inputString); % convert to ascii
    fn(fn<48)=[]; % remove non-alphanumeric characters
    fn=char(fn); % convert from uint8's to char
    outputStrings{stringIndex}=fn; % and store
end

if NStrings==1 % just one string passed?
    outputStrings=char(outputStrings); % then return char
end

outputStrings=genvarname(outputStrings);

end
