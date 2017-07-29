function [dfdt,dfldt,dfgdt,liquidBins,gasBins,internalVolumeBins,aggregationKernel,collisionEfficiency,collisionFrequency,breakageKernel]=Compartment(fComingIn,fgComingIn,fPreviousCompartment,flPreviousCompartment,fgPreviousCompartment,particleAverageVelocity,distanceBetweenCompartments,premixingTime,liquidAdditionTime,liquidAdditionRate,fl,fg,f,time,timeStep,aggregationKernelConstant,numberOfCollisions,breakageKernelConstant,consolidationConstant,minimumPorosity,granuleSaturationFactor,ns,nss,vs,vss,s_meshxy,ss_meshxy,saggregationCheck,ssaggregationCheck,sind,ssind,sind_b,ssind_b,s_high,s_low,ss_high,ss_low,s_check_b,ss_check_b,s_mesh_break,ss_mesh_break,nubmerOfDEMBins,timeStepDEM,collisionEfficiencyConstant,DEMDiameter,diameter,numberOfImpacts,breakageProbability)
%% INITIALIZATION
liquidBins(1:ns,1:nss)=0;
gasBins(1:ns,1:nss)=0;
internalLiquid(1:ns,1:nss)=0;
externalLiquid(1:ns,1:nss)=0;
externalLiquidContent(1:ns,1:nss)=0;
aggregationRate(1:ns,1:nss,1:ns,1:nss)=0;
depletionThroughAggregation(1:ns,1:nss)=0;
depletionOfGasThroughAggregation(1:ns,1:nss)=0;
depletionOfLiquidThroughAggregation(1:ns,1:nss)=0;
liquidBirthThroughAggregation(1:ns,1:nss)=0;
gasBirthThroughAggregation(1:ns,1:nss)=0;
firstSolidVolumeThroughAggregation(1:ns,1:nss)=0;
secondSolidVolumeThroughAggregation(1:ns,1:nss)=0;
birthThroughAggregation(1:ns,1:nss)=0;
firstSolidBirthThroughAggregation(1:ns,1:nss)=0;
secondSolidBirthThroughAggregation(1:ns,1:nss)=0;
breakageRate(1:ns,1:nss,1:ns,1:nss)=0;
depletionThroughBreakage(1:ns,1:nss)=0;
depletionOfGasThroughBreakage(1:ns,1:nss)=0;
depletionOfLiquidthroughBreakage(1:ns,1:nss)=0;
birthThroughBreakage1(1:ns,1:nss)=0;
birthThroughBreakage2(1:ns,1:nss)=0;
firstSolidBirthThroughBreakage(1:ns,1:nss)=0;
secondSolidBirthThroughBreakage(1:ns,1:nss)=0;
liquidBirthThroughBreakage2(1:ns,1:nss)=0;
gasBirthThroughBreakage2(1:ns,1:nss)=0;
liquidBirthThroughBreakage1(1:ns,1:nss)=0;
gasBirthThroughBreakage1(1:ns,1:nss)=0;
formationThroughBreakageCA(1:ns,1:nss)=0;
formationOfLiquidThroughBreakageCA(1:ns,1:nss)=0;
formationOfGasThroughBreakageCA(1:ns,1:nss)=0;
dfdt(1:ns,1:nss)=0;
transferThroughLiquidAddition(1:ns,1:nss)=0;
dfldt(1:ns,1:nss)=0;
transferThroughConsolidation(1:ns,1:nss)=0;
dfgdt(1:ns,1:nss)=0;
birth_agg_low_low(1:ns,1:nss)=0;
birth_agg_high_high(1:ns,1:nss)=0;
birth_agg_low_high(1:ns,1:nss)=0;
birth_agg_high_low(1:ns,1:nss)=0;
birth_agg_low_low_liq(1:ns,1:nss)=0;
birth_agg_high_high_liq(1:ns,1:nss)=0;
birth_agg_low_high_liq(1:ns,1:nss)=0;
birth_agg_high_low_liq(1:ns,1:nss)=0;
birth_agg_low_low_gas(1:ns,1:nss)=0;
birth_agg_high_high_gas(1:ns,1:nss)=0;
birth_agg_low_high_gas(1:ns,1:nss)=0;
birth_agg_high_low_gas(1:ns,1:nss)=0;

%% BIN
% Calculation of liquid and gas bins
for s=1:ns
    for ss=1:nss
        %         if and(f(s,ss)~=0,f(s,ss)>=1) % If there are particles in the bin...
        if f(s,ss)~=0 % If there are particles in the bin...
            liquidBins(s,ss)=fl(s,ss)/f(s,ss); % New liquid bins are calculated as total amount liquid in that size class divided by the number of particle in that size class
            gasBins(s,ss)=fg(s,ss)/f(s,ss); % New gas bins are calculated as total amount liquid in that size class divided by the number of particle in that size class
        else
            liquidBins(s,ss)=0; % No particles in the bin -> particles in that size class has no liquid
            gasBins(s,ss)=0; % No particles in the bin -> particles in that size class has no gas
        end
    end
end
% fg
% f
% gasBins
%Internal and external liquid demarcation
for s=1:ns
    for ss=1:nss
        internalLiquid(s,ss)=min(granuleSaturationFactor*gasBins(s,ss),liquidBins(s,ss));
        externalLiquid(s,ss)=max(0,liquidBins(s,ss)-internalLiquid(s,ss));
        externalLiquidContent(s,ss)=externalLiquid(s,ss)/liquidBins(s,ss);
    end
end
internalVolumeBins=s_meshxy+ss_meshxy+internalLiquid+gasBins;
externalVolumeBins=s_meshxy+ss_meshxy+liquidBins+gasBins;%%%%%%%%%%%%%%%%%%%%%%%
volumeBins=s_meshxy+ss_meshxy;%%%%%%%%%%%%%%%%%%%%%%%

%% AGGREGATION
% aggregationKernel=AggregationKernel(ns,nss,aggregationKernelConstant,aggregationKernelConstantALPHA,aggregationKernelConstantDELTA,aggregationKernelConstantGAMMA,externalVolumeBins,externalLiquid);
% aggregationKernel=DEMDependentAggregationKernel(aggregationKernelConstant,numberOfCollisions,f,timeStepDEM,collisionEfficiencyConstant,criticialExternalLiquid,externalLiquid,binderViscosity,coefficientOfRestitution,surfaceAsperity,averageHeightOfSurfaceLiquid,averageVelocityCompartment,meanParticleVelocityCompartment,standardDeviationParticleVelocityCompartment,ns,nss,nubmerOfDEMBins,DEMDiameter,diameter,timeStep)

[aggregationKernel,collisionEfficiency,collisionFrequency]=DEMDependentAggregationKernel(aggregationKernelConstant,numberOfCollisions,f,timeStepDEM,collisionEfficiencyConstant,ns,nss,nubmerOfDEMBins,DEMDiameter,diameter,timeStep,externalLiquidContent);

for s1=1:ns
    for ss1=1:nss
        for s2=1:ns
            for ss2=1:nss
                aggregationRate(s1,ss1,s2,ss2)=saggregationCheck(s1,ss1,s2,ss2)*ssaggregationCheck(s1,ss1,s2,ss2)*aggregationKernel(s1,ss1,s2,ss2)*f(s1,ss1)*f(s2,ss2);
            end
        end
    end
end
% aggregationRate
depletionThroughAggregation(1:ns,1:nss)=reshape(sum(sum(aggregationRate,1),2),[ns nss])+reshape(sum(sum(aggregationRate,3),4),[ns nss]);
for s=1:ns
    for ss=1:nss
        depletionOfGasThroughAggregation(s,ss)=depletionThroughAggregation(s,ss)*gasBins(s,ss);
        depletionOfLiquidThroughAggregation(s,ss)=depletionThroughAggregation(s,ss)*liquidBins(s,ss);
    end
end


% birthThroughAggregation=accumarray({sind(:),ssind(:)},aggregationRate(:));
% birthThroughAggregation=birthThroughAggregation(1:ns,1:nss);
%
% firstSolidBirthThroughAggregationFunction(1:ns,1:nss,1:ns,1:nss)=(a_mesh_a+a_mesh_b).*aggregationRate;
% firstSolidBirthThroughAggregation=accumarray({sind(:),ssind(:)},firstSolidBirthThroughAggregationFunction(:));
% firstSolidBirthThroughAggregation=firstSolidBirthThroughAggregation(1:ns,1:nss);
%
% secondSolidBirthThroughAggregationFunction(1:ns,1:nss,1:ns,1:nss)=(c_mesh_a+c_mesh_b).*aggregationRate;
% secondSolidBirthThroughAggregation=accumarray({sind(:),ssind(:)},secondSolidBirthThroughAggregationFunction(:));
% secondSolidBirthThroughAggregation=secondSolidBirthThroughAggregation(1:ns,1:nss);


for s1=1:ns
    for ss1=1:nss
        for s2=1:ns
            for ss2=1:nss
                for a=1:ns
                    for b=1:nss
                        if sind(s1,ss1,s2,ss2)==a
                            if ssind(s1,ss1,s2,ss2)==b
                                birthThroughAggregation(a,b)=birthThroughAggregation(a,b)+aggregationRate(s1,ss1,s2,ss2);
                                firstSolidBirthThroughAggregation(a,b)=firstSolidBirthThroughAggregation(a,b)+(vs(s1)+vs(s2))*aggregationRate(s1,ss1,s2,ss2);
                                secondSolidBirthThroughAggregation(a,b)=secondSolidBirthThroughAggregation(a,b)+(vss(ss1)+vss(ss2))*aggregationRate(s1,ss1,s2,ss2);
                                liquidBirthThroughAggregation(a,b)=liquidBirthThroughAggregation(a,b)+(liquidBins(s1,ss1)+liquidBins(s2,ss2))*aggregationRate(s1,ss1,s2,ss2);
                                gasBirthThroughAggregation(a,b)=gasBirthThroughAggregation(a,b)+(gasBins(s1,ss1)+gasBins(s2,ss2))*aggregationRate(s1,ss1,s2,ss2);
                            end
                        end
                    end
                end
            end
        end
    end
end

for s=1:ns
    for ss=1:nss
        if birthThroughAggregation(s,ss)~=0
            firstSolidVolumeThroughAggregation(s,ss)=firstSolidBirthThroughAggregation(s,ss)/birthThroughAggregation(s,ss);
            secondSolidVolumeThroughAggregation(s,ss)=secondSolidBirthThroughAggregation(s,ss)/birthThroughAggregation(s,ss);
        else
            firstSolidVolumeThroughAggregation(s,ss)=0;
            secondSolidVolumeThroughAggregation(s,ss)=0;
        end
    end
end
arep=repmat(vs',[1 nss]); %'
crep=repmat(vss,[ns 1]);


birth_agg_low_low(1:ns-1,1:nss-1)=(arep(2:ns,1:nss-1)-firstSolidVolumeThroughAggregation(1:ns-1,1:nss-1))./(arep(2:ns,1:nss-1)-arep(1:ns-1,1:nss-1)).*(crep(1:ns-1,2:nss)-secondSolidVolumeThroughAggregation(1:ns-1,1:nss-1))./(crep(1:ns-1,2:nss)-crep(1:ns-1,1:nss-1)).*birthThroughAggregation(1:ns-1,1:nss-1);
birth_agg_high_high(2:ns,2:nss)=(firstSolidVolumeThroughAggregation(1:ns-1,1:nss-1)-arep(1:ns-1,1:nss-1))./(arep(2:ns,1:nss-1)-arep(1:ns-1,1:nss-1)).*(secondSolidVolumeThroughAggregation(1:ns-1,1:nss-1)-crep(1:ns-1,1:nss-1))./(crep(1:ns-1,2:nss)-crep(1:ns-1,1:nss-1)).*birthThroughAggregation(1:ns-1,1:nss-1);
birth_agg_low_high(1:ns-1,2:nss)=(arep(2:ns,1:nss-1)-firstSolidVolumeThroughAggregation(1:ns-1,1:nss-1))./(arep(2:ns,1:nss-1)-arep(1:ns-1,1:nss-1)).*(secondSolidVolumeThroughAggregation(1:ns-1,1:nss-1)-crep(1:ns-1,1:nss-1))./(crep(1:ns-1,2:nss)-crep(1:ns-1,1:nss-1)).*birthThroughAggregation(1:ns-1,1:nss-1);
birth_agg_high_low(2:ns,1:nss-1)=(firstSolidVolumeThroughAggregation(1:ns-1,1:nss-1)-arep(1:ns-1,1:nss-1))./(arep(2:ns,1:nss-1)-arep(1:ns-1,1:nss-1)).*(crep(1:ns-1,2:nss)-secondSolidVolumeThroughAggregation(1:ns-1,1:nss-1))./(crep(1:ns-1,2:nss)-crep(1:ns-1,1:nss-1)).*birthThroughAggregation(1:ns-1,1:nss-1);

birth_agg_low_low_liq(1:ns-1,1:nss-1)=(arep(2:ns,1:nss-1)-firstSolidVolumeThroughAggregation(1:ns-1,1:nss-1))./(arep(2:ns,1:nss-1)-arep(1:ns-1,1:nss-1)).*(crep(1:ns-1,2:nss)-secondSolidVolumeThroughAggregation(1:ns-1,1:nss-1))./(crep(1:ns-1,2:nss)-crep(1:ns-1,1:nss-1)).*liquidBirthThroughAggregation(1:ns-1,1:nss-1);
birth_agg_high_high_liq(2:ns,2:nss)=(firstSolidVolumeThroughAggregation(1:ns-1,1:nss-1)-arep(1:ns-1,1:nss-1))./(arep(2:ns,1:nss-1)-arep(1:ns-1,1:nss-1)).*(secondSolidVolumeThroughAggregation(1:ns-1,1:nss-1)-crep(1:ns-1,1:nss-1))./(crep(1:ns-1,2:nss)-crep(1:ns-1,1:nss-1)).*liquidBirthThroughAggregation(1:ns-1,1:nss-1);
birth_agg_low_high_liq(1:ns-1,2:nss)=(arep(2:ns,1:nss-1)-firstSolidVolumeThroughAggregation(1:ns-1,1:nss-1))./(arep(2:ns,1:nss-1)-arep(1:ns-1,1:nss-1)).*(secondSolidVolumeThroughAggregation(1:ns-1,1:nss-1)-crep(1:ns-1,1:nss-1))./(crep(1:ns-1,2:nss)-crep(1:ns-1,1:nss-1)).*liquidBirthThroughAggregation(1:ns-1,1:nss-1);
birth_agg_high_low_liq(2:ns,1:nss-1)=(firstSolidVolumeThroughAggregation(1:ns-1,1:nss-1)-arep(1:ns-1,1:nss-1))./(arep(2:ns,1:nss-1)-arep(1:ns-1,1:nss-1)).*(crep(1:ns-1,2:nss)-secondSolidVolumeThroughAggregation(1:ns-1,1:nss-1))./(crep(1:ns-1,2:nss)-crep(1:ns-1,1:nss-1)).*liquidBirthThroughAggregation(1:ns-1,1:nss-1);

birth_agg_low_low_gas(1:ns-1,1:nss-1)=(arep(2:ns,1:nss-1)-firstSolidVolumeThroughAggregation(1:ns-1,1:nss-1))./(arep(2:ns,1:nss-1)-arep(1:ns-1,1:nss-1)).*(crep(1:ns-1,2:nss)-secondSolidVolumeThroughAggregation(1:ns-1,1:nss-1))./(crep(1:ns-1,2:nss)-crep(1:ns-1,1:nss-1)).*gasBirthThroughAggregation(1:ns-1,1:nss-1);
birth_agg_high_high_gas(2:ns,2:nss)=(firstSolidVolumeThroughAggregation(1:ns-1,1:nss-1)-arep(1:ns-1,1:nss-1))./(arep(2:ns,1:nss-1)-arep(1:ns-1,1:nss-1)).*(secondSolidVolumeThroughAggregation(1:ns-1,1:nss-1)-crep(1:ns-1,1:nss-1))./(crep(1:ns-1,2:nss)-crep(1:ns-1,1:nss-1)).*gasBirthThroughAggregation(1:ns-1,1:nss-1);
birth_agg_low_high_gas(1:ns-1,2:nss)=(arep(2:ns,1:nss-1)-firstSolidVolumeThroughAggregation(1:ns-1,1:nss-1))./(arep(2:ns,1:nss-1)-arep(1:ns-1,1:nss-1)).*(secondSolidVolumeThroughAggregation(1:ns-1,1:nss-1)-crep(1:ns-1,1:nss-1))./(crep(1:ns-1,2:nss)-crep(1:ns-1,1:nss-1)).*gasBirthThroughAggregation(1:ns-1,1:nss-1);
birth_agg_high_low_gas(2:ns,1:nss-1)=(firstSolidVolumeThroughAggregation(1:ns-1,1:nss-1)-arep(1:ns-1,1:nss-1))./(arep(2:ns,1:nss-1)-arep(1:ns-1,1:nss-1)).*(crep(1:ns-1,2:nss)-secondSolidVolumeThroughAggregation(1:ns-1,1:nss-1))./(crep(1:ns-1,2:nss)-crep(1:ns-1,1:nss-1)).*gasBirthThroughAggregation(1:ns-1,1:nss-1);

formationThroughAggregationCA=birth_agg_high_high+birth_agg_low_high+birth_agg_high_low+birth_agg_low_low;
formationOfLiquidThroughAggregationCA=birth_agg_high_high_liq+birth_agg_low_high_liq+birth_agg_high_low_liq+birth_agg_low_low_liq;
formationOfGasThroughAggregationCA=birth_agg_high_high_gas+birth_agg_low_high_gas+birth_agg_high_low_gas+birth_agg_low_low_gas;

%% BREAKAGE
% breakageKernel=BreakageKernel(ns,nss,breakageKernelConstant,shearRate,externalVolumeBins);

breakageKernel=DEMDependentBreakageKernel(breakageKernelConstant,ns,nss,f,numberOfImpacts,nubmerOfDEMBins,DEMDiameter,diameter,timeStep,breakageProbability,timeStepDEM);

for s1=1:ns
    for ss1=1:nss
        for s2=1:ns
            for ss2=1:nss
                breakageRate(s1,ss1,s2,ss2)= s_check_b(s1,ss1,s2,ss2)*ss_check_b(s1,ss1,s2,ss2)*breakageKernel(s1,ss1,s2,ss2)*f(s1,ss1);
                depletionThroughBreakage(s1,ss1)=depletionThroughBreakage(s1,ss1)+breakageRate(s1,ss1,s2,ss2);
                depletionOfLiquidthroughBreakage(s1,ss1)=depletionThroughBreakage(s1,ss1)*liquidBins(s1,ss1);
                depletionOfGasThroughBreakage(s1,ss1)=depletionThroughBreakage(s1,ss1)*gasBins(s1,ss1);
            end
        end
    end
end

for s1=1:ns
    for ss1=1:nss
        for s2=1:ns
            for ss2=1:nss
                birthThroughBreakage1(s1,ss1)=birthThroughBreakage1(s1,ss1)+breakageRate(s1,ss1,s2,ss2);
            end
        end
    end
end

for s1=1:ns
    for ss1=1:nss
        for s2=1:ns
            for ss2=1:nss
                for a=1:ns-1
                    for b=1:nss-1
                        if sind_b(s1,ss1,s2,ss2)==a
                            if ssind_b(s1,ss1,s2,ss2)==b
                                birthThroughBreakage2(a,b)=birthThroughBreakage2(a,b)+breakageRate(s1,ss1,s2,ss2);
                                firstSolidBirthThroughBreakage(a,b)=firstSolidBirthThroughBreakage(a,b)+s_mesh_break(s1,ss1,s2,ss2)*breakageRate(s1,ss1,s2,ss2);
                                secondSolidBirthThroughBreakage(a,b)=secondSolidBirthThroughBreakage(a,b)+ss_mesh_break(s1,ss1,s2,ss2)*breakageRate(s1,ss1,s2,ss2);
                                liquidBirthThroughBreakage2(a,b)=liquidBirthThroughBreakage2(a,b)+(liquidBins(s1,ss1)*(1-(volumeBins(s2,ss2)/volumeBins(s1,ss1))))*breakageRate(s1,ss1,s2,ss2);
                                gasBirthThroughBreakage2(a,b)=gasBirthThroughBreakage2(a,b)+(gasBins(s1,ss1)*(1-(volumeBins(s2,ss2)/volumeBins(s1,ss1))))*breakageRate(s1,ss1,s2,ss2);
                            end
                        end
                    end
                end
                liquidBirthThroughBreakage1(s2,ss2)=liquidBirthThroughBreakage1(s2,ss2)+(liquidBins(s1,ss1)*(volumeBins(s2,ss2)/volumeBins(s1,ss1)))*breakageRate(s1,ss1,s2,ss2);
                gasBirthThroughBreakage1(s2,ss2)=gasBirthThroughBreakage1(s2,ss2)+(gasBins(s1,ss1)*(volumeBins(s2,ss2)/volumeBins(s1,ss1)))*breakageRate(s1,ss1,s2,ss2);
            end
        end
    end
end

firstSolidVoumeThroughBreakage=firstSolidBirthThroughBreakage./birthThroughBreakage2;
firstSolidVoumeThroughBreakage(firstSolidVoumeThroughBreakage~=firstSolidVoumeThroughBreakage)=0;
secondSolidVolumeThroughBreakage=secondSolidBirthThroughBreakage./birthThroughBreakage2;
secondSolidVolumeThroughBreakage(secondSolidVolumeThroughBreakage~=secondSolidVolumeThroughBreakage)=0;

fractionBreakage00=(s_high-s_low-abs(s_low-firstSolidVoumeThroughBreakage(1:ns,1:nss)))./(s_high-s_low).*(ss_high-ss_low-abs(ss_low-secondSolidVolumeThroughBreakage(1:ns,1:nss)))./(ss_high-ss_low);
fractionBreakage01=(s_high-s_low-abs(s_low-firstSolidVoumeThroughBreakage(1:ns,1:nss)))./(s_high-s_low).*(ss_high-ss_low-abs(ss_high-secondSolidVolumeThroughBreakage(1:ns,1:nss)))./(ss_high-ss_low);
fractionBreakage10=(s_high-s_low-abs(s_high-firstSolidVoumeThroughBreakage(1:ns,1:nss)))./(s_high-s_low).*(ss_high-ss_low-abs(ss_low-secondSolidVolumeThroughBreakage(1:ns,1:nss)))./(ss_high-ss_low);
fractionBreakage11=(s_high-s_low-abs(s_high-firstSolidVoumeThroughBreakage(1:ns,1:nss)))./(s_high-s_low).*(ss_high-ss_low-abs(ss_high-secondSolidVolumeThroughBreakage(1:ns,1:nss)))./(ss_high-ss_low);

formationThroughBreakageCA(1:ns-1,1:nss-1)=formationThroughBreakageCA(1:ns-1,1:nss-1)+birthThroughBreakage2(1:ns-1,1:nss-1).*fractionBreakage00(1:ns-1,1:nss-1);
formationThroughBreakageCA(1:ns-1,2:nss)=formationThroughBreakageCA(1:ns-1,2:nss)+birthThroughBreakage2(1:ns-1,1:nss-1).*fractionBreakage01(1:ns-1,1:nss-1);
formationThroughBreakageCA(2:ns,1:nss-1)=formationThroughBreakageCA(2:ns,1:nss-1)+birthThroughBreakage2(1:ns-1,1:nss-1).*fractionBreakage10(1:ns-1,1:nss-1);
formationThroughBreakageCA(2:ns,2:nss)=formationThroughBreakageCA(2:ns,2:nss)+birthThroughBreakage2(1:ns-1,1:nss-1).*fractionBreakage11(1:ns-1,1:nss-1);

formationOfLiquidThroughBreakageCA(1:ns-1,1:nss-1)=formationOfLiquidThroughBreakageCA(1:ns-1,1:nss-1)+liquidBirthThroughBreakage2(1:ns-1,1:nss-1).*fractionBreakage00(1:ns-1,1:nss-1);
formationOfLiquidThroughBreakageCA(1:ns-1,2:nss)=formationOfLiquidThroughBreakageCA(1:ns-1,2:nss)+liquidBirthThroughBreakage2(1:ns-1,1:nss-1).*fractionBreakage01(1:ns-1,1:nss-1);
formationOfLiquidThroughBreakageCA(2:ns,1:nss-1)=formationOfLiquidThroughBreakageCA(2:ns,1:nss-1)+liquidBirthThroughBreakage2(1:ns-1,1:nss-1).*fractionBreakage10(1:ns-1,1:nss-1);
formationOfLiquidThroughBreakageCA(2:ns,2:nss)=formationOfLiquidThroughBreakageCA(2:ns,2:nss)+liquidBirthThroughBreakage2(1:ns-1,1:nss-1).*fractionBreakage11(1:ns-1,1:nss-1);

formationOfGasThroughBreakageCA(1:ns-1,1:nss-1)=formationOfGasThroughBreakageCA(1:ns-1,1:nss-1)+gasBirthThroughBreakage2(1:ns-1,1:nss-1).*fractionBreakage00(1:ns-1,1:nss-1);
formationOfGasThroughBreakageCA(1:ns-1,2:nss)=formationOfGasThroughBreakageCA(1:ns-1,2:nss)+gasBirthThroughBreakage2(1:ns-1,1:nss-1).*fractionBreakage01(1:ns-1,1:nss-1);
formationOfGasThroughBreakageCA(2:ns,1:nss-1)=formationOfGasThroughBreakageCA(2:ns,1:nss-1)+gasBirthThroughBreakage2(1:ns-1,1:nss-1).*fractionBreakage10(1:ns-1,1:nss-1);
formationOfGasThroughBreakageCA(2:ns,2:nss)=formationOfGasThroughBreakageCA(2:ns,2:nss)+gasBirthThroughBreakage2(1:ns-1,1:nss-1).*fractionBreakage11(1:ns-1,1:nss-1);

%%  Particle Transfer
particleMovement=fComingIn+fPreviousCompartment*((particleAverageVelocity*timeStep)/distanceBetweenCompartments).*(1-((s_meshxy+ss_meshxy)/max(max(s_meshxy+ss_meshxy))))-f*((particleAverageVelocity*timeStep)/distanceBetweenCompartments);
liquidMovement=flPreviousCompartment*((particleAverageVelocity*timeStep)/distanceBetweenCompartments).*(1-((s_meshxy+ss_meshxy)/max(max(s_meshxy+ss_meshxy))))-fl*((particleAverageVelocity*timeStep)/distanceBetweenCompartments);
gasMovement=fgComingIn+fgPreviousCompartment*((particleAverageVelocity*timeStep)/distanceBetweenCompartments).*(1-((s_meshxy+ss_meshxy)/max(max(s_meshxy+ss_meshxy))))-fg*((particleAverageVelocity*timeStep)/distanceBetweenCompartments);

% particleMovement=fComingIn+fPreviousCompartment*((particleAverageVelocity*timeStep)/distanceBetweenCompartments).*(1-((s_meshxy+ss_meshxy)/max(max(s_meshxy+ss_meshxy))))-f*((particleAverageVelocity*timeStep)/distanceBetweenCompartments).*(1-((s_meshxy+ss_meshxy)/max(max(s_meshxy+ss_meshxy))));
% liquidMovement=flPreviousCompartment*((particleAverageVelocity*timeStep)/distanceBetweenCompartments).*(1-((s_meshxy+ss_meshxy)/max(max(s_meshxy+ss_meshxy))))-fl*((particleAverageVelocity*timeStep)/distanceBetweenCompartments).*(1-((s_meshxy+ss_meshxy)/max(max(s_meshxy+ss_meshxy))));
% gasMovement=fgComingIn+fgPreviousCompartment*((particleAverageVelocity*timeStep)/distanceBetweenCompartments).*(1-((s_meshxy+ss_meshxy)/max(max(s_meshxy+ss_meshxy))))-fg*((particleAverageVelocity*timeStep)/distanceBetweenCompartments).*(1-((s_meshxy+ss_meshxy)/max(max(s_meshxy+ss_meshxy))));

%% RATE
if time>=premixingTime && time<=(premixingTime+liquidAdditionTime);
    liquidAdditionRate=liquidAdditionRate*timeStep; % Liquid spray rate in m3/sec
else
    liquidAdditionRate=0;
end

totalSolidVolume=0;
for s=1:ns
    for ss=1:nss
        totalSolidVolume=totalSolidVolume+f(s,ss)*(vs(s)+vss(ss));
    end
end

for s=1:ns
    for ss=1:nss
        dfdt(s,ss)=particleMovement(s,ss)+formationThroughAggregationCA(s,ss)-depletionThroughAggregation(s,ss)+birthThroughBreakage1(s,ss)+formationThroughBreakageCA(s,ss)-depletionThroughBreakage(s,ss);
        
        if totalSolidVolume~=0
            transferThroughLiquidAddition(s,ss)=liquidAdditionRate*((vs(s)+vss(ss))/totalSolidVolume);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            transferThroughLiquidAddition(s,ss)=0;
        end
        
%         dfldt(s,ss)=formationOfLiquidThroughAggregationCA(s,ss)-depletionOfLiquidThroughAggregation(s,ss)+liquidBirthThroughBreakage1(s,ss)+formationOfLiquidThroughBreakageCA(s,ss)-depletionOfLiquidthroughBreakage(s,ss);
        
        dfldt(s,ss)=liquidMovement(s,ss)+f(s,ss)*transferThroughLiquidAddition(s,ss)+formationOfLiquidThroughAggregationCA(s,ss)-depletionOfLiquidThroughAggregation(s,ss)+liquidBirthThroughBreakage1(s,ss)+formationOfLiquidThroughBreakageCA(s,ss)-depletionOfLiquidthroughBreakage(s,ss);
        if fg(s,ss)~=0
            %             transferThroughConsolidation(s,ss)=consolidationConstant*(rpm/1500)^rpmConsolidationConstant*(externalVolumeBins(s,ss)-externalLiquid(s,ss))*(1-minimumPorosity)/(vs(s)+vss(ss))*(gasBins(s,ss)-(minimumPorosity/(1-minimumPorosity))*(vs(s)+vss(ss))+internalLiquid(s,ss));
            transferThroughConsolidation(s,ss)=consolidationConstant*(internalVolumeBins(s,ss))*(1-minimumPorosity)/(vs(s)+vss(ss))*(gasBins(s,ss)-minimumPorosity/(1-minimumPorosity)*(vs(s)+vss(ss))+internalLiquid(s,ss));
        else
            transferThroughConsolidation(s,ss)=0;
        end
        dfgdt(s,ss)=gasMovement(s,ss)+f(s,ss)*transferThroughConsolidation(s,ss)+formationOfGasThroughAggregationCA(s,ss)-depletionOfGasThroughAggregation(s,ss)+gasBirthThroughBreakage1(s,ss)+formationOfGasThroughBreakageCA(s,ss)-depletionOfGasThroughBreakage(s,ss);
    end
end

return

