classdef PolyMesh
    % Class for processing and storing OpenFOAM polyMesh data
    
    % Anthony Ricciardi
    % January 2022
    
    properties
        boundary % [nBoundaries, 1 Boundary] boundary data
        faces % {nFaces,[nPointsPerFace uint32]} faces point indices (base 0)
        neighbour % [nNeighbours, 1 uint32] neighbour cell indices (base 0)
        owner % [nOwners, 1 uint32] owner cell indices (base 0)
        points % [nPoints,3 double] point locations
    end
    
    methods (Static = true)
        function obj = read()
            obj = PolyMesh;
            obj.boundary = Boundary.read();
            obj.faces = obj.readFaces();
            obj.neighbour = obj.readNeighbour();
            obj.owner = obj.readOwner();
            obj.points = obj.readPoints();
        end
    end
    methods (Static = true, Access = private)
        function faces = readFaces()
            % function to read polyMesh faces file
            fileString = fileread(fullfile('polyMesh','faces'));
            faceData = regexp(fileString,'\(\s*(?<points>[^\n()]+)\s*\)','names');
            nFaces = size(faceData,2);
            for i = nFaces:-1:1
                faces{i,1} = uint32( sscanf(faceData(i).points,'%d').' );
            end
        end
        function neighbour = readNeighbour()
            % function to read polyMesh neighbour file
            fileString = fileread(fullfile('polyMesh','neighbour'));
            elementList = regexp(fileString,'\((?<element>[^()])+\)','names');
            neighbour= uint32( sscanf(elementList.element,'%d') );
        end
        function owner = readOwner()
            % function to read polyMesh owner file
            fileString = fileread(fullfile('polyMesh','owner'));
            elementList = regexp(fileString,'\((?<element>[^()])+\)','names');
            owner = uint32( sscanf(elementList.element,'%d') );
        end
        function points = readPoints()
            fileString = fileread(fullfile('polyMesh','points'));
            pointData = regexp(fileString,'\s*\(\s*(?<x>[^ \s()]+)\s+(?<y>[^ \s()]+)\s+(?<z>[^ \s()]+)\s*\)','names');
            nPoints = size(pointData,2);
            points = zeros(nPoints,3);
            for i = 1:nPoints
                points(i,:) = [sscanf(pointData(i).x,'%f'),...
                              sscanf(pointData(i).y,'%f'),...
                              sscanf(pointData(i).z,'%f')];
            end
        end
    end
end

