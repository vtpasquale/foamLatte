% Class for isoparametric quadrilateral plate elements for thermal analysis
% 
% Anthony Ricciardi
%
classdef TeaCquad4 < Element
    
    properties
        eid % [uint32] Element identification number.
        pid % [uint32] Property identification number of a PBEAM entry.
        g   % [3,1 uint32] Grid point identification numbers of connection points [G1,G2,G3].
        
        % theta or mcid
        tFlag % [logical] specifies how t is used to define thickness of element
        t     % [4,1 double] Thickness of element at grid points G1 through G4 if TFLAG=0 or []. If TFLAG=1, thickness becomes a product of Ti and the thickness on the PSHELL card.
        
        gdof % [1,18 uint32] indices of global degrees of freedom associated with the element.
        R_eg % [18 x 18 double] rotation matrix from the element reference frame to the nodal displacement reference frame.
        k_e % [18 x 18 double] element stiffness matrix in the element reference frame
        m_e % [18 x 18 double] element mass matrix in the element reference frame
        
        area % [double] element area
        volume % [double] element volume
        mass % [double] element mass
        
% %         E2dMem % [3,3 double] Membrane elasticity matrix
% %         E2dBend % [3,3 double] Bending elasticity matrix
% %         E2dShear % [2,2 double] Transverse shear elasticity matrix
        
        % tNodes % [4,1 double] thickness of element at grid points G1 through G4
        T_e0 % [3,3 double] Transformation from basic to element reference frame
        x2D_e % [4,2 double] node xy locations in element reference frame
        
% %         % stress recovery data
% %         centerT % [double] element thickness at center point
% %         centerBm % [double] membrane strain-dispacement matrix at center point
% %         centerBp % [double] plate strain-dispacement matrix at center point
% %         
% %         isMembrane % [logical] true if property includes membrane terms
% %         isPlate % [logical] true if property includes plate terms
    end
    properties (Constant = true, Hidden = true)
        ELEMENT_TYPE = uint8(33); % [uint8] Element code corresponding to Nastran item codes documentation.
        VTK_CELL_CLASSNAME = 'VtkCellQuad'; % [char] Vtk cell classname
        VTK_CELL_TYPE = 9; % [uint8] VTK cell type number
        HDF5_ELEMENT_FORCE_CLASSNAME = 'Hdf5ElementForceQuad4';
        HDF5_STRAIN_CLASSNAME = 'Hdf5ElementStrainQuad4';
        HDF5_STRESS_CLASSNAME = 'Hdf5ElementStressQuad4';
        GAUSS_POINT = 1/sqrt(3);
        TEMP_DOF = uint8(1:6:24);
%         MEMBRANE_DOF = uint8([1,2,7,8,13,14,19,20]);
%         DRILLING_DOF = uint8([6,12,18,24]);
%         PLATE_DOF = uint8([3,4,5,9,10,11,15,16,17,21,22,23]);
        T_PSI_THETA = [...
     1     0     0     0     0     0     0     0     0     0     0     0
     0     0    -1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0    -1     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0    -1     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0    -1
     0     0     0     0     0     0     0     0     0     0     1     0];
        TZ_2_ELEMENT_DOF = [0 0 0 0
                            0 0 0 0
                            1 0 0 0
                           0 0 0 0
                           0 0 0 0
                           0 0 0 0
                            0 0 0 0
                            0 0 0 0
                            0 1 0 0
                           0 0 0 0
                           0 0 0 0
                           0 0 0 0
                            0 0 0 0
                            0 0 0 0
                            0 0 1 0
                           0 0 0 0
                           0 0 0 0
                           0 0 0 0
                            0 0 0 0
                            0 0 0 0
                            0 0 0 1
                           0 0 0 0
                           0 0 0 0
                           0 0 0 0];
    end
    methods
        function obj = TeaCquad4(cquad4)
            % Construct from Cquad4 object array
            if nargin > 0
                if ~isa(cquad4,'Cquad4'); error('Input must be type Cquad4'); end
                nElements = size(cquad4,1);
                for i = 1:nElements
                    teaCquad4 = TeaCquad4;
                    teaCquad4.eid = cquad4(i).eid;
                    teaCquad4.pid = cquad4(i).pid;
                    teaCquad4.g = cquad4(i).g;
                    teaCquad4.tFlag = cquad4(i).tFlag;
                    obj(i,1) = teaCquad4;
                end
            end
        end
        function obj=assemble_sub(obj,model)
            % geometry Data
            n1 = model.point.getNode(obj.g(1),model);
            n2 = model.point.getNode(obj.g(2),model);
            n3 = model.point.getNode(obj.g(3),model);
            n4 = model.point.getNode(obj.g(4),model);
            obj.gdof = [n1.gdof,n2.gdof,n3.gdof,n4.gdof];
            x1 = n1.x_0;
            x2 = n2.x_0;
            x3 = n3.x_0;
            x4 = n4.x_0;
            
            % element coordinate system
            x0 = .25*(x1+x2+x3+x4);
            xe = .5*(x2 + x3) - x0; xe = xe./normCS(xe);
            ze = cross3(x2-x0,x3-x0); ze = ze./normCS(ze);
            ye = cross3(ze, xe); ye = ye./normCS(ye);
            T_e0 = [xe,ye,ze].';
            
            % node positions in element coordinate system
            x_e = T_e0*([x1, x2, x3, x4] - [x0, x0, x0, x0]);
            x2D_e = x_e(1:2,:);
            
            % Property and material
            pshell = model.property.getProperty(obj.pid,model,'Pshell');
            pshellT = pshell.t;
% %             if pshell.isMembrane
% %                 obj.isMembrane = true;
% %                 obj.E2dMem = pshell.E2dMembrane;
% %             else
% %                 obj.isMembrane = false;
% %             end
% %             if pshell.isPlate
% %                 obj.isPlate = true;
% %                 obj.E2dBend = pshell.E2dBend;
% %                 obj.E2dShear = pshell.E2dShear;
% %             else
% %                 obj.isPlate = false;
% %             end


            % Hard code thermal conductivity
            k = 1;

            
            % Element thickness
            if obj.tFlag
                if isempty(pshellT); error('PSHELL thinkness must be specified if PSHELL TFLAG=0.'); end
                if isempty(obj.t); error('CQUAD4 thinkness must be specified if PSHELL TFLAG=0.'); end
                tNodes = pshellT*obj.t;
            else
                if ~isempty(obj.t)
                    if isempty(pshellT); error('No shell thinkness defined. Shell thickness must be specifed on PSHELL or CQUAD4 entries.'); end
                    tNodes = obj.t;
                else
                    tNodes = ones(4,1)*pshellT;
                end
            end
            
            % Four-point integration
            conductivityMatrix = zeros(4);
% %             if obj.isMembrane
% %                 membraneTranslationStiffness = zeros(8);
% %             end
% %             if obj.isPlate
% %                 bendingStiffness = zeros(12);
% %             end

            consistentMass = zeros(4);
            obj.volume = 0;
            obj.area = 0;
            Xi = obj.GAUSS_POINT*[-1  1  1 -1];
            Eta= obj.GAUSS_POINT*[-1 -1  1  1];
            for i = 1:4
                
                % Gauss point evaluation
                [N,dNdxi,dNdeta] = TeaCquad4.evaluateShapeFunctions(Xi(i),Eta(i));
                [J,detJ,invJ] =    TeaCquad4.calculateJacobian2D(dNdxi,dNdeta,x2D_e);
                
                % Thickness at point
                tGauss = N * tNodes;
                
                % Consistent mass matrix integration
                conductivityMatrix = conductivityMatrix + k*(dNdxi.'*dNdxi + dNdeta.'*dNdeta)*detJ;
                
                % Conductivity matrix integration
                consistentMass = consistentMass + N.'*N*(pshell.nsm+pshell.rho*tGauss)*detJ;
                
% %                 % Stiffness integration
% %                 if obj.isMembrane
% %                     Bm = Cquad4.calculateMembraneB(dNdxi,dNdeta,invJ);
% %                     membraneTranslationStiffness = membraneTranslationStiffness ...
% %                         + Bm(1:2,:).'*obj.E2dMem(1:2,1:2)*Bm(1:2,:)*detJ*tGauss;
% %                 end
% %                 if obj.isPlate
% %                     Bp = Cquad4.calculateMindlinPlateB(N,dNdxi,dNdeta,invJ);
% % %                     bendingStiffness = bendingStiffness ...
% % %                         + (tGauss^3/12)*Bp(1:2,:).'*obj.E2dBend(1:2,1:2)*Bp(1:2,:)*detJ;
% %                     bendingStiffness = bendingStiffness ...
% %                         + (tGauss^3/12)*Bp(1:3,:).'*obj.E2dBend*Bp(1:3,:)*detJ;
% %                 end

                % Volume and area integratino
                obj.volume = obj.volume + detJ*tGauss;
                obj.area = obj.area + detJ;
            end
            
% %             % Reduced-order integration and center point recovery data
% %             [N,dNdxi,dNdeta] = Cquad4.evaluateShapeFunctions(0,0);
% %             [J,detJ,invJ] =    Cquad4.calculateJacobian2D(dNdxi,dNdeta,x2D_e);
% %             obj.centerT = N * tNodes;
% %             if obj.isMembrane
% %                 obj.centerBm = Cquad4.calculateMembraneB(dNdxi,dNdeta,invJ);
% %                 membraneShearStiffness = 4*obj.centerBm(3,:).'*obj.E2dMem(3,3)*obj.centerBm(3,:)*detJ*obj.centerT;
% %             end
% %             if obj.isPlate
% %                 obj.centerBp = Cquad4.calculateMindlinPlateB(N,dNdxi,dNdeta,invJ);
% %                 transverseShearStiffness = 4*obj.centerBp(4:5,:).'*obj.E2dShear*obj.centerBp(4:5,:)*detJ*obj.centerT;
% % %                 bendingStiffness = bendingStiffness + 4*(tGauss^3/12)*obj.centerBp(3,:).'*obj.E2dBend(3,3)*obj.centerBp(3,:)*detJ;
% %             end
            
            % element total mass
            obj.mass = sum(consistentMass(:));
            
% %             % element stiffness matrix
% %             obj.k_e = zeros(24);
% %             if obj.isMembrane
% %                 obj.k_e(obj.MEMBRANE_DOF,obj.MEMBRANE_DOF) = membraneTranslationStiffness + membraneShearStiffness;
% % 
% %                 % drilling stiffness
% %                 % Nastran Element User's guide equation for CQUAD8 and CTRIA6
% %                 % krot = model.k6rot/1e6 * (obj.E2dBend(1,1)+obj.E2dBend(2,2))
% %                 
% %                 % reverse engineered equation for CQUAD4 K(6,6)
% %                 krot = model.k6rot/1e7 * obj.E2dMem(3,3);
% %                 
% %                 % Application similar to Zienkiewicz approach as documented
% %                 % in 16.4-2 in CMPW Eq. 16.4-2.
% %                 obj.k_e(obj.DRILLING_DOF,obj.DRILLING_DOF) = krot*(1/3)*[ 3,-1,-1,-1;
% %                                                                          -1, 3,-1,-1;
% %                                                                          -1,-1, 3,-1;
% %                                                                          -1,-1,-1, 3];
% %             end
% %             if obj.isPlate
% %                 obj.k_e(obj.PLATE_DOF,obj.PLATE_DOF) = bendingStiffness + transverseShearStiffness;
% %             end          
            
            % element conductivity matrix
            obj.k_e = zeros(24);
            obj.k_e(obj.TEMP_DOF,obj.TEMP_DOF) = conductivityMatrix;

            % element mass matrix
            me = zeros(24);
            me(1:6:24,1:6:24) = consistentMass;
            me(2:6:24,2:6:24) = consistentMass;
            me(3:6:24,3:6:24) = consistentMass;
            % lumped mass matrix
            if ~model.coupledMassFlag
                me = diag(sum(me,2));
            end
            obj.m_e = me;
            
            % Transformation matrix
            obj.R_eg(22:24,22:24) = T_e0*n4.T_g0.';
            obj.R_eg(19:21,19:21) = T_e0*n4.T_g0.';
            obj.R_eg(16:18,16:18) = T_e0*n3.T_g0.';
            obj.R_eg(13:15,13:15) = T_e0*n3.T_g0.';
            obj.R_eg(10:12,10:12) = T_e0*n2.T_g0.';
            obj.R_eg(7:9,7:9)     = T_e0*n2.T_g0.';
            obj.R_eg(4:6,4:6)     = T_e0*n1.T_g0.';
            obj.R_eg(1:3,1:3)     = T_e0*n1.T_g0.';
            
            % save select properties
            obj.T_e0 = T_e0;
            obj.x2D_e = x2D_e;
        end
        function [gdof,p_g]=processPressureLoad_sub(obj,pload)
            Xi = obj.GAUSS_POINT*[-1  1  1 -1];
            Eta= obj.GAUSS_POINT*[-1 -1  1  1];
            pressureForce_e = zeros(4,1);
            for i = 1:4 
                % Gauss point evaluation
                [N,dNdxi,dNdeta] = Cquad4.evaluateShapeFunctions(Xi(i),Eta(i));
                [~,detJ] =    Cquad4.calculateJacobian2D(dNdxi,dNdeta,obj.x2D_e);
                                
                % pressure at point
                pGauss = N.' .* pload.p;
                
                % Pressure integration
                pressureForce_e = pressureForce_e + detJ*pGauss;
            end
            p_e = obj.TZ_2_ELEMENT_DOF*pressureForce_e;
            p_g = obj.R_eg.'*p_e;
            gdof = obj.gdof; 
        end
        function [force,stress,strain,strainEnergy,kineticEnergy] = recover_sub(obj,u_g,model,returnFlags)
            % INPUTS
            % u_g [nGodf,nVectors double] Response vector in nodal displacement reference frame
            % returnFlags [1,5 logical] [force,stress,strain,strainEnergy,kineticEnergy] 1 -> recover, 0 -> return empty array []
            %
            % OUTPUTS
            % force = [14,nVectors double] Element forces
            %   indices:
            %    [1 |  Membrane force x
            %     2 |  Membrane force y
            %     3 |  Membrane force xy
            %     4 |  Bending moment x
            %     5 |  Bending moment y
            %     6 |  Bending moment xy
            %     7 |  Shear x
            %     8 |  Shear y           ]
            %
            % stress  = [8,nVectors double] Element stresses
            % strain  = [8,nVectors double] Element strains
            %   indices:
            %    [1 |  Z1=Fiber Distance 1
            %     2 |  Normal x at Z1
            %     3 |  Normal y at Z1
            %     4 |  Shear xy at Z1
            %     5 |  Shear angle at Z1
            %     6 |  Major principal at Z1
            %     7 |  Minor principal at Z1
            %     8 |  von Mises or maximum shear at Z1
            %     9 |  Z2=Fiber Distance 2
            %    10 |  Normal x at Z2
            %    11 |  Normal y at Z2
            %    12 |  Shear xy at Z2
            %    13 |  Shear angle at Z2
            %    14 |  Major principal at Z2
            %    15 |  Minor principal at Z2
            %    16 |  von Mises or maximum shear at Z2    ]
            %
            % strainEnergy  = [3,nVectors double] Element strain energy
            % kineticEnergy = [3,nVectors double] Element kinetic energy
            %   indices:
            %    [ energy
            %      energy----------> converted to percent total later by Element.recover()
            %      energyDensity];
            %  kineticEnergy scaled by omega later by Element.recover()
            
            % Check inputs
            if ~any(returnFlags); error('This function is not intended to be called if no vaules are to be recovered'); end
            
            % Element displacements
            u_e = obj.R_eg*u_g(obj.gdof,:);
            nVectors = size(u_e,2);
            
            % Membrane values
            if obj.isMembrane
                membraneStrain = obj.centerBm*u_e(obj.MEMBRANE_DOF,:);
                membraneForce = obj.centerT*obj.E2dMem*membraneStrain;
            else
                membraneStrain = zeros(3,nVectors);
                membraneForce = zeros(3,nVectors);
            end
            
            % Plate values
            if obj.isPlate
                plateCurvature = obj.centerBp*u_e(obj.PLATE_DOF,:);
                plateMoment = (obj.centerT^3/12)*obj.E2dBend*plateCurvature(1:3,:);
                plateShear = obj.E2dShear*plateCurvature(4:5,:);
            else
                plateCurvature = zeros(5,nVectors);
                plateMoment = zeros(3,nVectors);
                plateShear = zeros(2,nVectors);
            end
            
            % Force
            if returnFlags(1)
                force = zeros(8,nVectors);
                force(1:3,:) = membraneForce;
                force(4:6,:) = plateMoment;
                force(7:8,:) = plateShear;
            else
                force = [];
            end
            
            % stress or strain data
            if any(returnFlags(2:3))
                z = obj.centerT/2;
            end
            
            % Stress
            if returnFlags(2)
                membraneStress = (1/obj.centerT)*membraneForce;
                bendingStress = (6/obj.centerT^2)*plateMoment; 
                topStress    = membraneStress - bendingStress;
                bottomStress = membraneStress + bendingStress;
                
                vonMisesT = calculateVonMises(topStress);
                vonMisesB = calculateVonMises(bottomStress);
                
                [s1T,s2T,angleT] = calculatePrincipal(topStress);
                [s1B,s2B,angleB] = calculatePrincipal(bottomStress);
                
                stress = [
                    -z*ones(1,nVectors);
                    bottomStress(1,:)
                    bottomStress(2,:)
                    bottomStress(3,:)
                    angleB
                    s1B
                    s2B
                    vonMisesB
                    z*ones(1,nVectors);
                    topStress(1,:)
                    topStress(2,:)
                    topStress(3,:)
                    angleT
                    s1T
                    s2T
                    vonMisesT];
            else
                stress = [];
            end
            
            % Strain
            if returnFlags(3)
                topStrain    = membraneStrain - z*plateCurvature(1:3,:);
                bottomStrain = membraneStrain + z*plateCurvature(1:3,:);
                
                vonMisesT = calculateVonMises(topStrain);
                vonMisesB = calculateVonMises(bottomStrain);
                
                [s1T,s2T,angleT] = calculatePrincipal(topStrain);
                [s1B,s2B,angleB] = calculatePrincipal(bottomStrain);
                
                strain = [
                    -z*ones(1,nVectors);
                    bottomStrain(1,:)
                    bottomStrain(2,:)
                    bottomStrain(3,:)
                    angleB
                    s1B
                    s2B
                    vonMisesB
                    z*ones(1,nVectors);
                    topStrain(1,:)
                    topStrain(2,:)
                    topStrain(3,:)
                    angleT
                    s1T
                    s2T
                    vonMisesT];
            else
                strain = [];
            end
            
            % Strain Energy
            if returnFlags(4)
                strainEnergy0 = .5*diag(u_e.'*obj.k_e*u_e).';
                strainEnergy = [strainEnergy0;
                    strainEnergy0;%---> converted to percent total later by Element.recover()
                    (1/obj.volume)*strainEnergy0];
            else
                strainEnergy = [];
            end
            
            % Kinetic Energy (scaled by omega later by Element.recover())
            if returnFlags(5)
                kineticEnergy0 = .5*diag(u_e.'*obj.m_e*u_e).';
                kineticEnergy = [kineticEnergy0;
                    kineticEnergy0;%---> converted to percent total later by Element.recover()
                    (1/obj.volume)*kineticEnergy0];
            else
                kineticEnergy = [];
            end
        end
    end
    methods (Access=private, Static=true)
        function [N,dNdxi,dNdeta] = evaluateShapeFunctions(xi,eta)
            N      = .25*[(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)];
            dNdxi  = .25*[      -(1-eta),       (1-eta),       (1+eta),      -(1+eta)];
            dNdeta = .25*[(1-xi)*-1     ,(1+xi)*-1     ,(1+xi)        ,(1-xi)        ];
        end
        function [J,detJ,invJ]=calculateJacobian2D(dNdxi,dNdeta,x2D_e)
            J = [dNdxi;dNdeta]*x2D_e.';
            detJ = det(J);
            invJ=(1/detJ)*[J(2,2),-J(1,2);-J(2,1),J(1,1)];
        end
    end
end