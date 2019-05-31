% Compile SDependenceCalculator and mexInterface
%
% Files needed:
% mexSDependenceCalculator.cpp    % mex interface
% SDependenceCalculator.cpp       % Core algorythm
% SDependenceCalculator.h
% ObjectHandle.h                  % a handle class for mex-functions written by Tim Bailey (2004)

mex  -g mexSDependenceCalculator.cpp SDependenceCalculator.cpp
%% Calculate several quantities and the dependency of the probabilities on S

% Specify terminal segment numbers
NRange=1:23; % N is bound by memory availability, we found that we can go upto N=23 with 4 GB of internal memory
NSize=length(NRange);

% Specify S values
SStep= 0.25;
SRange=-1:SStep:1;

SDCHandle= mexSDependenceCalculator();
% Set maximum number of terminal segments
tic
mexSDependenceCalculator(SDCHandle, 'N',uint32(NRange(end)));
toc
PCell=cell(3,1);
for N=1:NSize
    treeCount=mexSDependenceCalculator(SDCHandle, 'C',uint32(NRange(N)));
    PCell{N}=zeros(treeCount,length(SRange));
end



for S=1:length(SRange)
    % Set and evaluate S values
    mexSDependenceCalculator(SDCHandle, 'S',double(SRange(S)))
    
    % Put S-dependent probabilities in cell object PCell
    for N=1:NSize        
        % Calculate probabilities    
        PCell{N}(:,S)=mexSDependenceCalculator(SDCHandle, 'P',uint32(NRange(N)));
    end
end

for N=1:NSize        
        filename=['mat-files/SDCRes_N',num2str(NRange(N)),'_SRange',num2str(100*SRange(1)),'_',num2str(100*SStep),'_',num2str(100*SRange(end))]
        A=mexSDependenceCalculator(SDCHandle, 'A',uint32(NRange(N)));
        C=mexSDependenceCalculator(SDCHandle, 'C',uint32(NRange(N)));
        E=mexSDependenceCalculator(SDCHandle, 'E',uint32(NRange(N)));
        H=mexSDependenceCalculator(SDCHandle, 'H',uint32(NRange(N)));
        M=mexSDependenceCalculator(SDCHandle, 'M',uint32(NRange(N)));
        P=PCell{N};
        save(filename,'N','C','A','E','H','M','P','SStep','SRange')
end
    
mexSDependenceCalculator(SDCHandle, 'D' ); % Free memory
