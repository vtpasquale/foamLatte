classdef VolumeField
    % A field defined at cell centers and boundary face centeroids
    
    properties
        name % [char]
        value % [nCells+nBoundaryFaces,dimension double] field values. See MMD Fig. 7.5 for array structure.
    end
    methods
        function obj = VolumeField(nameIn,valueIn)
            if nargin > 0
                obj.name = nameIn;
                obj.value = valueIn;
            end
        end
        function faceValue = interpolateValueToInteriorFaces(obj,fvMesh)
            
            % get the list of owners and neighbours for all interior faces
            interiorFaces = 1:fvMesh.boundary(1).startFace;
            interiorOwners = fvMesh.owner(interiorFaces)+1;
            interiorNeighbours = fvMesh.neighbour(interiorFaces)+1;
            
            % get the geometric factor (interpolation factor)for all interior faces
            gf = fvMesh.interpFactor(interiorFaces);
                        
            % compute the linear interpolation of phi into all interior faces
            faceValue =   obj.value(interiorNeighbours,:) .* gf ...
                        + obj.value(interiorOwners,:) .*  (1-gf);
        end
        function grad = computeGaussGradient(obj,fvMesh)
            % Compute the Gauss Gradient
            % Following Moukalled, Mangani, and Darwish (2015)
            % [nValues,nDim]=size(obj.value);
            nCells=size(fvMesh.cellFaces,1);
            
            % Initialize gradValue Array
            gradValue = zeros(nCells,3);
            
            %% INTERIOR FACES contribution to gradient
            % get the list of interior faces indices and the list of boundary faces indices.
            firstBoundaryFaceIndex = fvMesh.boundary(1).startFace + 1;
            lastBoundaryFaceIndex = fvMesh.boundary(end).startFace + fvMesh.boundary(end).nFaces;
            interiorFaces = 1:(firstBoundaryFaceIndex-1);
            boundaryFaces = firstBoundaryFaceIndex:lastBoundaryFaceIndex;
            nBoundaryFaces = lastBoundaryFaceIndex-firstBoundaryFaceIndex+1;
            clear firstBoundaryFaceIndex lastBoundaryFaceIndex
                       
            % get the list of element indices and the list of boundary element indices
            boundaryValueIndices = (nCells+1):(nCells+nBoundaryFaces);
            
            % get the list of owners and neighbours for all interior faces
            interiorOwners = fvMesh.owner(interiorFaces)+1;
            interiorNeighbours = fvMesh.neighbour(interiorFaces)+1;
            % get the surface vector of all interior faces
            Sf = fvMesh.Sf(interiorFaces,:);
                        
            % compute the linear interpolation of phi into all interior faces
            faceValue = obj.interpolateValueToInteriorFaces(fvMesh);
                     
            % loop over faces and add contribution of face flux to the owner and neighbour elements of the faces
            % tricky to vectorize
            for i=interiorFaces
                gradValue(interiorOwners(i),:) = gradValue(interiorOwners(i),:) + faceValue(i)*Sf(i,:);
                gradValue(interiorNeighbours(i),:) = gradValue(interiorNeighbours(i),:) - faceValue(i)*Sf(i,:);
            end
                        
            %% BOUNDARY FACES contribution to gradient
            % get the list of elements owning a boundary face
            boundaryFaceOwner = fvMesh.owner(boundaryFaces)+1;
            % get the boundary values
            boundaryValue = obj.value(boundaryValueIndices,:);
            % get the surface vector of all boundary faces
            Sb = fvMesh.Sf(boundaryFaces,:);
            % loop over all boundary faces and add phi flux contribution to gradient
            for k=1:nBoundaryFaces
                gradValue(boundaryFaceOwner(k),:) = gradValue(boundaryFaceOwner(k),:) + boundaryValue(k)*Sb(k,:);
            end
            
            %% Divide by cell volume
            gradValue = gradValue./fvMesh.cellVolume;
            
            %% Store in VolumeField 
            grad = VolumeField(['grad_',obj.name],gradValue);

        end
        function cellDataToVtk(obj,fvMesh,fid)
            % print cell data to VTK file
            nCells=size(fvMesh.cellFaces,1);
            switch size(obj.value,2)
                case 1
                    fprintf(fid,'SCALARS %s float 1\n',obj.name);
                    fprintf(fid,'LOOKUP_TABLE default\n');
                    fprintf(fid,'%f\n',obj.value(1:nCells));
                case 3
                    fprintf(fid,'VECTORS %s float\n',obj.name);
                    fprintf(fid,'%f %f %f\n',obj.value(1:nCells,:).');
                otherwise
                    error('Update function for field dimension.')
            end
        end
    end
    methods (Static=true)
        function obj = constructFromEquation(nameIn,fvMesh,functionHandle)
            % construct field from functionHandle = @(x,z,y)
            nDimensions = size( functionHandle(0,0,0), 2);
            nCells = size(fvMesh.cellFaces,1);
            nValues = fvMesh.boundary(end).startFace + fvMesh.boundary(end).nFaces;
            value = zeros(nValues,nDimensions);
            value(1:nCells,:) = functionHandle(fvMesh.cellCentroid(:,1),fvMesh.cellCentroid(:,2),fvMesh.cellCentroid(:,3));
            nBoundaries = size(fvMesh.boundary,1);
            for i = 1:nBoundaries
                faceIndices = (1:fvMesh.boundary(i).nFaces)+ fvMesh.boundary(i).startFace;
                dofIndices  = nCells + (faceIndices-fvMesh.boundary(1).startFace);
                value(dofIndices,:) = functionHandle(fvMesh.faceCentroid(faceIndices,1),fvMesh.faceCentroid(faceIndices,2),fvMesh.faceCentroid(faceIndices,3));
            end
            obj = VolumeField(nameIn,value);
        end
    end
end

