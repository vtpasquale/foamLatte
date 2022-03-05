clear all; close all; clc
addpath(genpath('C:\Users\vtpas\Documents\MATLAB\NASTRAN_CoFE\NASTRAN_CoFE\nastran_cofe\'));


cofe = Cofe('obliqueParallelogramCoFE.dat','assemble',false);

%% Process entities
% elements
elements = cofe.model.element;
nElements = size(elements,1); 
if ~isa(elements,'Cquad4'); error('Modify for more element types.'); end

% nodes
nodes = cofe.model.point;
nNodes = size(nodes,1); 
if ~isa(nodes,'Node'); error('Scalor points not allowed.'); end

%% Convert elements to thermal
% teaCquad4 = TeaCquad4(elements);
% cofe.model.element = teaCquad4;
fvCquad4 = FvCquad4(elements);
cofe.model.element = fvCquad4;

%% CoFE SOL 101 Process
% Assemble model
cofe.model = cofe.model.assemble();

% Static solution presolve
cofe.model(1).reducedModel=cofe.model(1).reducedModel.solveUaAllLoadSets();

% Solve and recover
cofe.solution = Solution.constructFromModel(cofe.model);
cofe.solution = cofe.solution.solve(cofe.model);

% % Results
% hdf5 = Hdf5(cofe);
% hdf5.export(['obliqueParallelogram','.h5']);

vtkFile = VtkFile(cofe);
vtkFile.print('fvObliqueParallelogram.vtk')