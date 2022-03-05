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
teaCquad4 = TeaCquad4(elements);
fvCquad4 = FvCquad4(elements);

fe = cofe;
fv = cofe;

fe.model.element = teaCquad4;
fv.model.element = fvCquad4;

fe.model = fe.model.assemble();
fv.model = fv.model.assemble();

fek = fe.model.element(200).k_e(1:6:end,1:6:end)
fvk = fv.model.element(200).k_e(1:6:end,1:6:end)

fek./fvk