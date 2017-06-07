clear all
close all
clc

extension = '*.csv'; % stores the extension of the file that needs to be searched for
nameoffilesfound = dir(extension); % the files which exist in the folder with this extension

names = {nameoffilesfound.name,1}; % storing names in an array

for i = 1:length(names)-1 % it includes the 1 at the end for some reason
    flag = 0;
    flag2 = 0;
    A1 = csvread(string(names(i)),1,0);
    A1diff = zeros(length(A1(:,1))-1,17);
    Asum = zeros(length(A1(:,1))-1,1);
    for j = 1:length(A1(:,1))-1
        for m = 1:17
            A1diff(j,m) = sum(abs(A1(i+1,m) - A1(i,m)));
        end
        Asum(j) = sum(A1(j,[2:end]));
    end
    for k = 1:10000:length(A1diff)
        %for l = 2:17
            if (sum(A1diff(k,[2:end]))/Asum(k) < 0.1 & sum(A1diff(i+1,[2:end]))/Asum(k) < 0.1 & sum(A1diff(k+2,[2:end]))/Asum(k) < 0.1 & sum(A1diff(k+3,[2:end]))/Asum(k) < 0.1 & sum(A1diff(k+4,[2:end]))/Asum(k) < 1)
                flag =1; k1 = k; 
            else
                flag = 0;
            end
            
        %end
    end
    
    %for l = 2:17
        for r = 1:1000:length(A1diff)-1000
            if (sum(A1diff(r,[2:end]))/Asum(r) < 0.1 & sum(A1diff(r+200,[2:end]))/Asum(r) < 0.1 & sum(A1diff(r+400,[2:end]))/Asum(r) < 0.1 & sum(A1diff(r+600,[2:end]))/Asum(r) < 0.1 & sum(A1diff(r+800,[2:end]))/Asum(r) < 1 & A1diff(r,2)~=0)
                flag2 =1; k2 = r; A1diff(r,:)
                break
            else
                flag2 = 0;
            end
        end
    %end
    
    
    if (flag == 1)
        name2 = cell2mat(names(i));
        disp(sprintf('%s is at steady state at Time = %d',name2([1:end-4]),A1(k1,1)));
    end
    
    if (flag2 == 1)
        name2 = cell2mat(names(i));
        disp(sprintf('%s is at steady state at Time = %d over a longer period of time',name2([1:end-4]),A1(k2,1)));
    end
    
end
