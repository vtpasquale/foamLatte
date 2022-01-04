clear all; close all; clc

fvMesh = FiniteVolumeMesh(PolyMesh.read());

fvMesh.writeVTK('test.vtk')
