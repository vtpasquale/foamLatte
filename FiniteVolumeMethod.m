classdef FiniteVolumeMethod
    
    % Anthony Ricciardi
    % February 2022
    
    properties
        gamma % diffusion coefficent
        fluxC1 % [nFaces,1 double] owner cell flux coefficent
        fluxC2 % [nFaces,1 double] neighbour cell flux coefficent
        fluxV % [nFaces,1 double] non-orthogonal flux term
        fluxT % [nFaces,1 double] total flux passing through face
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
            obj.fluxC1(interiorFaces,1) = obj.gamma .* gDiffInterior;
            obj.fluxC2(interiorFaces,1) =-obj.gamma .* gDiffInterior;
            obj.fluxV(interiorFaces,1) =  obj.gamma .* dot(gradInteriorFaceValue,TfInterior,2);
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
                TPatch = fvMesh.T(faceIndices,:);
                ownerPatch = fvMesh.owner(faceIndices)+1;
                obj.fluxC1(faceIndices) = obj.gamma.*gDiffPatch;
                obj.fluxC2(faceIndices) = - obj.gamma.*gDiffPatch;
                obj.fluxV(faceIndices) =  - obj.gamma.*dot(gradBoundaryPatch,TPatch,2);
            end
            
            
        end
        
    end
    
end

