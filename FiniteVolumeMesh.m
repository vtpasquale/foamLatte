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
    properties (Hidden = true) % from processFaceGeometry()
        faceCentroid % [nFaces,3 double] face centroids
        Sf % [nFaces,3 double] face normals
        faceArea % [nFaces,1 double] face areas
    end
    properties (Hidden = true) % from processCellGeometry()
        cellCentroid % [nCells,3 double] cell centroids
        cellVolume % [nCells,1 double] cell volumes
    end
    properties (Hidden = true) % from calculateFaceQuantities()
        interpFactor % [nInteriorFaces,1 double] face interpolation factor
        gDiff % [nFaces,1 double] geometric diffusion coefficient
        T % [nFaces,3 double] Nonorthogonal face normal components
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
        function obj = processFaceGeometry(obj)
            % Process face geometry
            % Following Moukalled, Mangani, and Darwish (2015)
            nFaces = size(obj.faces,1);            
            obj.faceCentroid = zeros(nFaces,3);
            obj.Sf = zeros(nFaces,3);
            obj.faceArea = zeros(nFaces,1);
            for i=1:nFaces
                facePoints = obj.faces{i} + 1;
                nFacePoints = size(facePoints,2);
                
                % a point is constructed within the polygon based on the
                % average of all the points that define the polygon.
                % This point is the geometric centre of the polygon xG
                geometricCenter = mean(obj.points(facePoints,:),1);
                
                % Using the geometric centre (Fig. 6.18) a number of 
                % triangles are formed with the centre as the apex for each
                % side of the polygon. For each of the triangles, the 
                % centroid (since for triangles xG and xCE coincide) and 
                % area are readily computed.
                centroid = [0, 0, 0];
                Sf_ = [0, 0, 0];
                area = 0;
                %
                % using the center compute the area and centroid
                % of virtual triangles based on the centre and the
                % face nodes
                %
                for j=1:nFacePoints
                    point1 = geometricCenter;
                    point2 = obj.points( facePoints(j) ,:);
                    if(j<nFacePoints)
                        point3 = obj.points( facePoints(j+1) ,:);
                    else
                        point3 = obj.points( facePoints(1) ,:);
                    end
                    localCentroid = (point1+point2+point3)/3;
                    localSf = 0.5*cross(point2-point1,point3-point1);
                    localArea = sqrt(localSf*localSf.');
                    centroid = centroid + localArea*localCentroid;
                    Sf_ = Sf_ + localSf;
                    area = area + localArea;
                end
                centroid = centroid/area;
                %
                obj.faceCentroid(i,:) = centroid;
                obj.Sf(i,:) = Sf_;
                obj.faceArea(i) = area;
            end
        end
        function obj = processCellGeometry(obj)
            % Process cell geometry
            % Following Moukalled, Mangani, and Darwish (2015)
            
            % The process starts by computing the location of
            % the geometric centre of the polyhedron element and 
            % decomposing it into a number of polygonal pyramids.           
            nCells = size(obj.cellFaces,1);
            obj.cellCentroid = zeros(nCells,3);
            obj.cellVolume = zeros(nCells,1);
            for i = 1:nCells
                
                % Compute a rough cell center
                center = mean(obj.points(obj.cellPoints{i}+1,:),1);
                
                % using the centre, compute the area and centroid
                % of virtual triangles based on the centre and the
                % face nodes
                cellFaces_ = obj.cellFaces{i}+1;
                nCellFaces = size(cellFaces_,2);
                localVolumeCentroidSum = [0, 0, 0];
                localVolumeSum = 0;
                for j=1:nCellFaces
                    Sf_ = obj.Sf(cellFaces_(j),:);
                    Cf = obj.faceCentroid(cellFaces_(j),:) - center;
                    % calculate face-pyramid volume
                    localVolume = Sf_*Cf.'/3;
                    if localVolume < 0
                        localVolume = -1* localVolume;
                    end
                    % Calculate face-pyramid centre
                    localCentroid = 0.75*obj.faceCentroid(cellFaces_(j),:) + 0.25*center;
                    %Accumulate volume-weighted face-pyramid centre
                    localVolumeCentroidSum = localVolumeCentroidSum + localCentroid*localVolume;
                    % Accumulate face-pyramid volume
                    localVolumeSum = localVolumeSum + localVolume;
                end
                centroid = localVolumeCentroidSum/localVolumeSum;
                volume = localVolumeSum;
                %
                obj.cellCentroid(i,:) = centroid;
                obj.cellVolume(i) = volume;
            end
            
        end
        function obj = calculateFaceQuantities(obj)
            % Calculate face interpolation factor (weighting factor)
            % Following Moukalled, Mangani, and Darwish (2015)
            nFaces = size(obj.faces,1);
            nInteriorFaces = size(obj.neighbour,1);
            nBoundaryFaces = nFaces - nInteriorFaces;
            obj.interpFactor = zeros(nInteriorFaces,1);
            
            % Interior face quantities
            obj.gDiff=zeros(nFaces,1);
            obj.T=zeros(nFaces,3);
            for i=1:nInteriorFaces
                ef = obj.Sf(i,:)./norm(obj.Sf(i,:));
                dCf = obj.faceCentroid(i,:)-obj.cellCentroid(obj.owner(i)+1,:);
                dfF = obj.cellCentroid(obj.neighbour(i)+1,:) - obj.faceCentroid(i,:);
                obj.interpFactor(i) = dCf*ef.' / (dCf*ef.' + dfF*ef.');

                dCF = obj.cellCentroid(obj.neighbour(i)+1,:)-obj.cellCentroid(obj.owner(i)+1,:);
                ndCF = norm(dCF);
                eCF = dCF./ndCF;
                E = obj.faceArea(i)*eCF;
                obj.gDiff(i) = norm(E)/ndCF;
                obj.T(i,:) = obj.Sf(i,:) - E;
            end
            
            % Bounary Face Quantities
            % This approach is deduced; not documented.
            for i=(nInteriorFaces+1):nBoundaryFaces
                dCf = obj.faceCentroid(i,:)-obj.cellCentroid(obj.owner(i)+1,:);
                ndCf = norm(dCf);
                eCF = dCf./ndCf;
                E = obj.faceArea(i)*eCF;
                obj.gDiff(i) = norm(E)/ndCF;
                obj.T(i,:) = obj.Sf(i,:) - E;
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
            fprintf(fid,'CELLS %d %d\n',nCells,sum(obj.nPointsPerCell)+nCells);
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

