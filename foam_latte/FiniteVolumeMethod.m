classdef FiniteVolumeMethod
    
    % Anthony Ricciardi
    % February 2022
    
    properties
        gamma % diffusion coefficent
        fluxC1 % [nFaces,1 double] owner cell flux coefficent
        fluxC2 % [nFaces,1 double] neighbour cell flux coefficent
        fluxV % [nFaces,1 double] non-orthogonal flux term
        fluxT % [nFaces,1 double] total flux passing through face
        
        A % [nCells,nCells sparse double] coefficent matrix
        b % [nCells,1      sparse double] RHS vector
    end
    
    methods
        function obj = FiniteVolumeMethod(fvMesh,gammaIn,volField)
            nFaces = size(fvMesh.faces,1);
            obj.gamma = gammaIn;
            obj.fluxC1 = zeros(nFaces,1);
            obj.fluxC2 = zeros(nFaces,1);
            obj.fluxV  = zeros(nFaces,1);
            obj.fluxT = zeros(nFaces,1);
            
            %% Process field
            gradVolField = volField.computeGaussGradient(fvMesh);
            gradInteriorFaceValue = gradVolField.interpolateValueToInteriorFaces(fvMesh);
            
            %% Assembly of the diffusion term at interior faces
            interiorFaces = 1:fvMesh.boundary(1).startFace;
            interiorOwner = fvMesh.owner(interiorFaces)+1;
            interiorNeighbour = fvMesh.neighbour(interiorFaces)+1;
            interiorFaces = 1:fvMesh.boundary(1).startFace;
            gDiffInterior = fvMesh.gDiff(interiorFaces);
            TfInterior = fvMesh.T(interiorFaces,:);
            obj.fluxC1(interiorFaces,1) =   obj.gamma .* gDiffInterior;
            obj.fluxC2(interiorFaces,1) = - obj.gamma .* gDiffInterior;
            obj.fluxV(interiorFaces,1) =  - obj.gamma .* dot(gradInteriorFaceValue,TfInterior,2);
            obj.fluxT(interiorFaces,1) = ...
                obj.fluxC1(interiorFaces) .* volField.value(interiorOwner) ...
              + obj.fluxC2(interiorFaces) .* volField.value(interiorNeighbour) ...
              + obj.fluxV(interiorFaces);
          
            %% Assembly of the diffusion term at boundary faces
            nBoundaries = size(fvMesh.boundary,1);
            for i = 1:nBoundaries
                boundaryPatch = fvMesh.boundary(i);
                faceIndices = boundaryPatch.startFace + (1:boundaryPatch.nFaces);
                gDiffPatch = fvMesh.gDiff(faceIndices);
                % TPatch = fvMesh.T(faceIndices,:);
                ownerPatch = fvMesh.owner(faceIndices)+1;
                obj.fluxC1(faceIndices) =   obj.gamma.*gDiffPatch;
                % obj.fluxC2(faceIndices) = - obj.gamma.*gDiffPatch;
                
                % specified volField value
                obj.fluxV(faceIndices)  = zeros(boundaryPatch.nFaces,1);
                % - obj.gamma.*dot(gradBoundaryPatch,TPatch,2);
                
                obj.fluxT(faceIndices,1) = ...
                      obj.fluxC1(faceIndices) .* volField.value(ownerPatch) ...
                    + obj.fluxC2(faceIndices) .* volField.value(faceIndices) ...
                    + obj.fluxV(faceIndices);
                
            end
            
            %% Assemble 
            nCells = size(fvMesh.cellFaces,1);
            a = SparseTriplet(3*nCells);
            b = SparseTriplet(nCells);
            
            % Fluxes of interior faces
            for i = interiorFaces
                iOwner = fvMesh.owner(i)+1;
                iNeighbour = fvMesh.neighbour(i)+1;
                gDof = [iOwner;iNeighbour];
                aLocal = [ obj.fluxC1(i)  obj.fluxC2(i);
                          -obj.fluxC1(i) -obj.fluxC2(i)];
                bLocal = [-obj.fluxT(i);
                           obj.fluxT(i)];
                a.addMatrix(aLocal,gDof);
                b.addVector(bLocal,gDof);
            end
            
            % Fluxes of boundary faces
            allBoundaryFaces = (fvMesh.boundary(1).startFace+1):nFaces;
            for i=allBoundaryFaces
                gDof = fvMesh.owner(i)+1;
                aLocal =   obj.fluxC1(i);
                bLocal = - obj.fluxT(i);
                a.addMatrix(aLocal,gDof);
                b.addVector(bLocal,gDof);
            end
            
            % Convert triplit
            obj.A = a.convertToSparseMatrix(nCells,nCells);
            obj.b = b.convertToSparseMatrix(nCells,1);
        end
        
    end
    
end

