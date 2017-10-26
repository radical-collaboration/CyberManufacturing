%% DEM Dependent Breakage Kernel
function breakageKernel=DEMDependentBreakageKernel(breakageKernelConstant,ns,nss,f,numberOfImpacts,nubmerOfDEMBins,DEMDiameter,diameter,timeStep,breakageProbability,timeStepDEM)

%% INITIALIZATION
impactFrequency(1:ns,1:nss)=0;
breakageKernel(1:ns,1:nss,1:ns,1:nss)=0;

%% Impact Frequency (from 1D Number of Impacts)
scaledDEMDiameter=DEMDiameter*(max(max(diameter))/max(max(DEMDiameter)));

for s1=1:ns
    for ss1=1:nss
        for i=1:nubmerOfDEMBins-1
            if f(s1,ss1)>=1
                if diameter(s1,ss1) <=scaledDEMDiameter(i+1) && diameter(s1,ss1)>scaledDEMDiameter(i)
                    impactFrequency(s1,ss1)=(numberOfImpacts(i)*timeStep)/(f(s1,ss1)*timeStepDEM);
                end
            end
        end
    end
end

%% Breakage Kernel Calculation
for s1=1:ns
    for ss1=1:nss
        for s2=1:ns
            for ss2=1:nss
                breakageKernel(s1,ss1,s2,ss2)=breakageKernelConstant*impactFrequency(s1,ss1)*breakageProbability;
            end
        end
    end
end

return

