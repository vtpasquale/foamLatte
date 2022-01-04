classdef FiniteVolumeMesh < PolyMesh
    % Class for processing and storing finite volume mesh data
    
    % Anthony Ricciardi
    % January 2022
    
    properties
        cellFaces % {nCells,[nFacesPerCell uint32]} cell face indices (base 0)
        cellPoints % {nCells,[nPointsPerCell uint32]} cell vertex point indices (base 0)
    end
    properties (Hidden = true)
        nFacesPerCell % [nCells,1 uint32] number of faces per cell
        nPointsPerCell % [nCells,1 uint32] number of faces per cell
    end
    
    methods
        function obj = FiniteVolumeMesh(polyMesh)
            if ~isa(polyMesh,'PolyMesh')
                error('The input should be type PolyMesh')
            end
            obj.boundary = polyMesh.boundary;
            obj.faces    = polyMesh.faces;
            obj.neighbour= polyMesh.neighbour;
            obj.owner    = polyMesh.owner;
            obj.points   = polyMesh.points;
            
            % create cell connectivity from polyMesh data
            nCells = max(obj.neighbour)+1;
            obj.cellFaces = cell(nCells,1);
            obj.cellPoints = cell(nCells,1);
            obj.nFacesPerCell = zeros(nCells,1,'uint32');
            obj.nPointsPerCell = zeros(nCells,1,'uint32');
            for i = uint32(1):uint32(nCells)
                obj.cellFaces{i} = [ find(obj.owner==(i-1))
                    find(obj.neighbour==(i-1)) ].' - 1;
                obj.nFacesPerCell(i) = size(obj.cellFaces{i},2);
                if size(unique(obj.cellFaces{i}))~=obj.nFacesPerCell(i)
                    error('There is a mesh issue')
                end
                
                
                if obj.nFacesPerCell(i) ~= uint32(6)
                    error('Update to handle more element types')
                else
                    obj.nPointsPerCell(i) = uint32(8);
                    
                    % Assemble local faces temporarily
                    localFace = zeros(6,4,'uint32');
                    localFace(1,:) =  obj.faces{ obj.cellFaces{i}(1) + 1 };
                    for j = 2:obj.nFacesPerCell(i)
                        nextFacePoints =  obj.faces{ obj.cellFaces{i}(j) + 1 };
                        if all( ismember(localFace(1,1:2),nextFacePoints) )
                            localFace(2,:) =  obj.faces{ obj.cellFaces{i}(j) + 1 };
                        elseif all( ismember(localFace(1,2:3),nextFacePoints) )
                            localFace(3,:) =  obj.faces{ obj.cellFaces{i}(j) + 1 };
                        elseif all( ismember(localFace(1,3:4),nextFacePoints) )
                            localFace(4,:) =  obj.faces{ obj.cellFaces{i}(j) + 1 };
                        elseif all( ismember(localFace(1,[4,1]),nextFacePoints) )
                            localFace(5,:) =  obj.faces{ obj.cellFaces{i}(j) + 1 };
                        else
                            if any( ismember(localFace(1,:),nextFacePoints) )
                                error('Mesh processing issue')
                            end
                            localFace(6,:) =  obj.faces{ obj.cellFaces{i}(j) + 1 };
                        end
                    end
                    
                    % order cell vertex points in VTK format order
                    vertex4 = localFace(2,...
                        ismember(localFace(2,:),localFace(5,:)) & ...
                        ismember(localFace(2,:),localFace(6,:))   );
                    vertex5 = localFace(2,...
                        ismember(localFace(2,:),localFace(3,:)) & ...
                        ismember(localFace(2,:),localFace(6,:))   );
                    vertex6 = localFace(4,...
                        ismember(localFace(4,:),localFace(3,:)) & ...
                        ismember(localFace(4,:),localFace(6,:))   );
                    vertex7 = localFace(4,...
                        ismember(localFace(4,:),localFace(5,:)) & ...
                        ismember(localFace(4,:),localFace(6,:))   );
                    obj.cellPoints{i} = [localFace(1,:),vertex4,vertex5,vertex6,vertex7];
                    if size(unique(obj.cellPoints{i}))~=size(obj.cellPoints{i})
                        error('There is a mesh issue')
                    end
                end
            end
        end
        function writeVTK(obj,filename)
            % Write the mesh to a .vtk file.
            % https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
            fid = fopen(filename,'w+');
            fprintf(fid,'# vtk DataFile Version 3.0\n');
            fprintf(fid,'My example\n');
            fprintf(fid,'ASCII\n');
            
            fprintf(fid,'\nDATASET UNSTRUCTURED_GRID\n');
            fprintf(fid,'POINTS %d float\n',size(obj.points,1));
            fprintf(fid,'%f %f %f\n',obj.points.');
            
            nCells = size(obj.cellFaces,1);
            fprintf(fid,'CELLS %d %d\n',nCells,nCells+sum(obj.nPointsPerCell) );
            for i = 1:nCells
                fprintf(fid,'%d ',obj.nPointsPerCell(i));
                fprintf(fid,'%d ',obj.cellPoints{i});
                fprintf(fid,'\n');
            end
            
            cellTypes = zeros(nCells,1);
            cellTypes(  obj.nFacesPerCell==uint32(6)) = uint32(12); % VolHexs
            if any(cellTypes == int32(0) )
                error('Update to handle other cell types')
            end
            fprintf(fid,'CELL_TYPES %d\n',nCells);
            fprintf(fid,'%d\n',cellTypes);            
            fclose(fid);
        end
    end
end

