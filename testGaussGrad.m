clear all; close all; clc

fvMesh = FiniteVolumeMesh(PolyMesh.read());
fvMesh = fvMesh.processFaceGeometry();
fvMesh = fvMesh.processCellGeometry();
fvMesh = fvMesh.calculateFaceInterpolationFactor();
fvMesh = fvMesh.calculateGeometricDiffusionCoefficient();

nCells = size(fvMesh.cellFaces,1);

functionHandle = @(x,y,z) 10*x.^2.*y.^2;
gradFunctionHandle = @(x,y,z) [20*x.* y.^2, 20*y.* x.^2, 0.0.*x];
T = VolumeField.constructFromEquation('T',fvMesh,functionHandle);
gradT = VolumeField.constructFromEquation('gradT',fvMesh,gradFunctionHandle);
gaussGradT = VolumeField('gaussGradT',T.computeGaussGradient(fvMesh));
gaussGradError = VolumeField('gaussGradError',gaussGradT.value-gradT.value(1:nCells,:));

fvMesh.writeVTK('test.vtk')
fid = fopen('test.vtk','a+');
fprintf(fid,'CELL_DATA %d\n',nCells);
T.cellDataToVtk(fvMesh,fid)
gradT.cellDataToVtk(fvMesh,fid)
gaussGradT.cellDataToVtk(fvMesh,fid)
gaussGradError.cellDataToVtk(fvMesh,fid)
fclose(fid);

