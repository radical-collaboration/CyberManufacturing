%% DEM Dependent Aggregation Kernel
function [aggregationKernel,collisionEfficiency,collisionFrequency]=DEMDependentAggregationKernel(aggregationKernelConstant,numberOfCollisions,f,timeStepDEM,collisionEfficiencyConstant,ns,nss,nubmerOfDEMBins,DEMDiameter,diameter,timeStep,externalLiquidContent)
% function aggregationKernel=DEMDependentAggregationKernel(aggregationKernelConstant,numberOfCollisions,f,timeStepDEM,collisionEfficiencyConstant,criticialExternalLiquid,externalLiquid,binderViscosity,coefficientOfRestitution,surfaceAsperity,averageHeightOfSurfaceLiquid,averageVelocityCompartment,meanParticleVelocityCompartment,standardDeviationParticleVelocityCompartment,ns,nss,nubmerOfDEMBins,DEMDiameter,diameter)

%% Initialization
collisionFrequency(1:ns,1:nss,1:ns,1:nss)=0;
collisionEfficiency(1:ns,1:nss,1:ns,1:nss)=0;
aggregationKernel(1:ns,1:nss,1:ns,1:nss)=0;

%% Collision Frequency
% for s1=1:ns
%     for ss1=1:nss
%         for s2=1:ns
%             for ss2=1:nss
%                 collisionFrequency(s1,ss1,s2,ss2)=numberOfCollisions(s1,ss1,s2,ss2)/(f(s1,ss1)*f(s2,ss2)*timeStepDEM);
%             end
%         end
%     end
% end

%% Collision Frequency (from 2D Number of Collisions)
scaledDEMDiameter=DEMDiameter*(max(max(diameter))/max(max(DEMDiameter)));
% mid=min(min(diameter))
% mad=max(max(diameter))
% midd=min(min(scaledDEMDiameter))
% madd=max(max(scaledDEMDiameter))

collisionFrequency(1:ns,1:nss,1:ns,1:nss)=1;
% for s1=1:ns
%     for ss1=1:nss
%         for s2=1:ns
%             for ss2=1:nss
%                 for i=1:nubmerOfDEMBins-1
%                     for j=1:nubmerOfDEMBins-1
%                         if f(s1,ss1)>=0 && f(s2,ss2)>=0
%                             if diameter(s1,ss1) <=scaledDEMDiameter(i+1) && diameter(s1,ss1)>scaledDEMDiameter(i)
%                                 if diameter(s2,ss2) <=scaledDEMDiameter(j+1) && diameter(s2,ss2)>scaledDEMDiameter(j)
%                                     collisionFrequency(s1,ss1,s2,ss2)=(numberOfCollisions(i,j)*timeStep)/(f(s1,ss1)*f(s2,ss2)*timeStepDEM);
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

% MCF=max(max(max(max(collisionFrequency))))

%% Constant Collision Efficiency
% FROM Barrasso, Tamrakar, Ramachandran. Procedia Engineering 102 (2015) 1295–1304. (p. 1298)
% collisionEfficiency(1:ns,1:nss,1:ns,1:nss)=collisionEfficiencyConstant;

%%  Constant Collision Efficiency Value Applied Based on External Liquid Content
% FROM Barrasso, Ramachandran. ChERD 93 (2015) 304–317. (p. 308)
% collisionEfficiencyConstant=0.001;
% criticialExternalLiquid=0.15;

% FROM Sen, Barrasso, Singh, Ramachandran. Processes 2014, 2, 89-111. (p. 96)
collisionEfficiencyConstant=0.01;
criticialExternalLiquid=0.2;

for s1=1:ns
    for ss1=1:nss
        for s2=1:ns
            for ss2=1:nss
                if externalLiquidContent(s1,ss1)>=criticialExternalLiquid && externalLiquidContent(s2,ss2)>=criticialExternalLiquid
                    collisionEfficiency(s1,ss1,s2,ss2)=collisionEfficiencyConstant;
                else
                    collisionEfficiency(s1,ss1,s2,ss2)=0;
                end
            end
        end
    end
end

%% Collision Efficiency Calculated from Particle Veloctiy Values
% FROM Barrasso, Ramachandran. J Pharm Innov 2016. (p. ___) AND Barrasso, Eppinger, Pereira, Aglave, Debus, Bermingham, Ramachandran. CES 123 (2015) 500–513. (p. 504)
% diameterHarmonicMean
% massHarmonicMean
% criticalVelocity=((3*pi*diameterHarmonicMean^2*binderViscosity)/(8*massHarmonicMean))*(1+1/coefficientOfRestitution)*log(averageHeightOfSurfaceLiquid/surfaceAsperity);

% muVelocityCompartment=functionOf(meanParticleVelocityCompartment,standardDeviationParticleVelocityCompartment);
% sigmaVelocityCompartment=functionOf(meanParticleVelocityCompartment,standardDeviationParticleVelocityCompartment);
% velocityDistributionCompartment=(1/(averageVelocityCompartment*sqrt(2*pi)*sigmaVelocityCompartment))*exp(-(((log(averageVelocityCompartment)-muVelocityCompartment)^2)/(2*sigmaVelocityCompartment^2)));

% collisionEfficiency(s1,ss1,s2,ss2)=functionof(criticalVelocity,velocityDistributionCompartment);

%% Aggregation Kernel Calculation
for s1=1:ns
    for ss1=1:nss
        for s2=1:ns
            for ss2=1:nss
                aggregationKernel(s1,ss1,s2,ss2)=aggregationKernelConstant*collisionFrequency(s1,ss1,s2,ss2)*collisionEfficiency(s1,ss1,s2,ss2);
            end
        end
    end
end
% MAK=max(max(max(max(aggregationKernel))))

return

