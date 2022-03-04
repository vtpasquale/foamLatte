classdef Boundary
    % Class for processing and storing OpenFOAM boundary data
    
    % Anthony Ricciardi
    % January 2022
    
    properties
        name % [char] boundary name
        type % [char] boundary mesh type
        nFaces % [uint32] number of faces
        startFace % [uint32] boundary start face (base 0)
    end
    properties (SetAccess = private)
        conditionType % [char] bounary condition type
    end
    methods (Static = true)
        function obj = read()
            fileString = fileread(fullfile('constant','polyMesh','boundary'));
            splitFile = regexp(fileString,'\((?<patches>[^()])+\)','names');
            
            splitNames = regexp(splitFile.patches,'(?<patchNames>[^{}])+\s*\{','names');
            splitKeyValues = regexp(splitFile.patches,'\{(?<keyValues>[^{}])+\}','names');
                        
            % Check interim results
            nNames = size(splitNames,2);
            nKeyValueBlocks = size(splitKeyValues,2);
            if nNames ~= nKeyValueBlocks
                error('Error parsing boundary file');
            end
            nBoundaries = nNames;
            
            % Initialize
            obj(nBoundaries,1) = Boundary;
            
            % Find the key values
            for i = 1:nBoundaries
                obj(i).name = strtrim( splitNames(i).patchNames );
                
                % type
                typeString = regexp(splitKeyValues(i).keyValues,...
                    'type\s*(?<type>\S)+\s*;','names');
                if ~all( size(typeString)==[1,1] )
                    error('Error parsing boundary type');
                end
                obj(i).type = typeString.type;
                
                % nFaces
                nFacesString = regexp(splitKeyValues(i).keyValues,...
                    'nFaces\s*(?<nFaces>\S)+\s*;','names');
                if ~all( size(nFacesString)==[1,1] )
                    error('Error parsing boundary type');
                end
                obj(i).nFaces = uint32( sscanf(nFacesString.nFaces,'%d') );
                               
                % startFace
                startFaceString = regexp(splitKeyValues(i).keyValues,...
                    'startFace\s*(?<startFace>\S)+\s*;','names');
                if ~all( size(startFaceString)==[1,1] )
                    error('Error parsing boundary type');
                end
                obj(i).startFace = uint32( sscanf(startFaceString.startFace,'%d') );
            end
        end
    end
    methods
        function obj = setConditionType(obj,name,conditionTypeIn)
            % Sets condition type for 'name' in Boundary array
            [n,m]=size(obj);
            if n < 1; error('Can''t set condition type for empty Boundary array.'); end
            if m ~=1; error('Boundary array m ~=1.'); end
            names = {obj.name};
            index = strcmp(names,name);
            if sum(index)~=1; error('One and only one boundary patch should match name: %s',name); end
            obj(index).conditionType = conditionTypeIn;
        end        
    end
    
end

