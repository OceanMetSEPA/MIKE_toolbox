function modelInfoStruct = modelRunSpecifications( mikeSpecificationFile )
% Create struct with parameters used in MIKE specification file (e.g. .m21fm)
%
% All settings used within model run are specified within struct (you just
% need to find them!)
%
% MIKE File consists of nested blocks with format e.g.
% [DOMAIN]
% ... various settings
% EndSect // DOMAIN
% These are nested accordingly within output struct
%
% INPUTS:
% mikeSpecificationFile (e.g. .m21fm, .m3fm)
%
% OUTPUT:
% struct containing filename and model settings
%
% This function determines the nesting of each block and generates a string
% to identify the nested structure where the data should be stored. The
% data is then assigned using 'eval'. There is probably a better way of
% doing this, but it seems to work...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   modelRunSpecifications.m  $
% $Revision:   1.0  $
% $Author:   ted.schlicke  $
% $Date:   Jul 21 2016 13:27:40  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin~=1
    help MIKE.modelRunSpecifications
    return
end

if ~exist(mikeSpecificationFile,'file')
    error('File ''%s'' not found',mikeSpecificationFile)
end

txt=readTxtFile(mikeSpecificationFile,'trim',1,'empty',0);

% First of all, find block sections. These are lines which start with '['
% after discarding white space.
i0Blocks=cellfun(@(x)x(1)=='[',txt);
i0Blocks=find(i0Blocks);
blockNames=strtrim(txt(i0Blocks));
blockNames=cellfun(@(x)x(2:end-1),blockNames,'Unif',0);
numberOfBlocks=length(blockNames);
% Find end of these blocks. For each block start index, we check the rest
% of the file for the corresponding 'EndSect' index:
i1Blocks=arrayfun(@(i)stringFinder(txt(i0Blocks(i):end),{'EndSect',blockNames{i}},'output','index','first',1),1:numberOfBlocks)'+i0Blocks-1;

if numberOfBlocks==0
    error('No blocks found in file- is it a MIKE model file?')
elseif length(i1Blocks)~=length(i0Blocks)
    error('Blocks not correctly identified')
end

% Store these in a struct
blockInfoStruct=struct2struct(struct('Blocks',blockNames,'StartIndex',num2cell(i0Blocks),'EndIndex',num2cell(i1Blocks)));

% Now determine nesting depth:
nestDepths=NaN(numberOfBlocks,1);
nestTxt=cell(numberOfBlocks,1);
for blockIndex=1:numberOfBlocks
    i0=i0Blocks(blockIndex);
    % Nesting depth: add one for every section start; subtract one for
    % every section end. Add 1 for good luck (actually so root level = 1)
    ithNestDepth=sum(i0Blocks<i0)-sum(i1Blocks<i0)+1; %
    nestDepths(blockIndex)=ithNestDepth;
    parentIndices=NaN(1,ithNestDepth);
    for nestIndex=ithNestDepth:-1:1
        parentIndex=find(nestDepths(1:blockIndex)==(nestIndex),1,'last');
        parentIndices(nestIndex)=parentIndex;
    end
    % Prepare a string with nested levels separated by '.'. 
    nestedField=strcat(blockNames(parentIndices),'.');
    nestedField=horzcat(nestedField{:});
    nestedField(end)=[]; % remove last '.'    
    nestTxt{blockIndex}=nestedField; 
end
% Add these to struct
blockInfoStruct.Nest=nestDepths;
blockInfoStruct.Path=nestTxt;

% Fill in structure. We do this by extracting text corresponding to each
% block range (removing nested blocks). These are converted to a struct of
% values. This struct is then added to our master struct using the nest
% string 

% Remove base nest (FemEngineHD)
blockInfoStruct.Path=strrep(blockInfoStruct.Path,strcat(blockInfoStruct.Path{1},'.'),'');

modelInfoStruct=struct('File',mikeSpecificationFile);
maxNest=max(nestDepths);
for nestIndex=2:maxNest % Start at 2 to ignore base nest
    nestStruct=structFilter(blockInfoStruct,blockInfoStruct.Nest==nestIndex);
    NNestStructs=length(nestStruct.Blocks);
    nestStruct=struct2struct(nestStruct);
    for structIndex=1:NNestStructs
        iNestStruct=nestStruct(structIndex) ;
        i0=iNestStruct.StartIndex;
        i1=iNestStruct.EndIndex;
        % Extract text that's between start & end indices, but not within
        % nested blocks:
        k=[i0:min(blockInfoStruct.StartIndex(blockInfoStruct.StartIndex>i0)),max(blockInfoStruct.EndIndex(blockInfoStruct.EndIndex<i1)):i1];
        itxt=txt(k);
        parameterStruct=splitStringByEqualsSign(itxt);
        if rand(1)>Inf,disp(parameterStruct),end % never true - but suppresses warning about unused 'parameterStruct' variable (which is actually used below)
        cmd=sprintf('modelInfoStruct.%s=parameterStruct',iNestStruct.Path);
        evalc(cmd);
    end
end
