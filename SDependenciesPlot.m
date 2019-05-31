%% Set parameter ranges used in calculations
SStep= 0.25;
SRange=-1:SStep:1;

N=11;

% load calculated results
filename=['mat-files/SDCRes_N',num2str(N),'_SRange',num2str(100*SRange(1)),'_',num2str(100*SStep),'_',num2str(100*SRange(end))];
load(filename,'A','C','E','H','M', 'P','SStep','SRange')

%Define bins for  summed electronic path lengths
Nbin=35  %#ok<NOPTS> % Number of bins for summed electronic pathlengths

edge1=min(E)-0.01*(max(E)-min(E));
edgeN=max(E)+0.01*(max(E)-min(E));
edges=edge1:(edgeN-edge1)/Nbin:edgeN;
centers=(edges(2:end)+edges(1:end-1))/2;
[n,bin] = histc(E,edges);  %#ok<ASGLU>

%Define bins for tree asymmetry indices
NAbin=20 %#ok<NOPTS> % Number of bins for tree asymmetrie index
Aedge1=min(A)-0.01*(max(A)-min(A));
AedgeN=max(A)+0.01*(max(A)-min(A));
Aedges=Aedge1:(AedgeN-Aedge1)/NAbin:AedgeN;
Acenters=(Aedges(2:end)+Aedges(1:end-1))/2;
[n,Abin] = histc(A,Aedges);

%% Set plot parameters
% markers and their sizes
markers={'ko-','ko-','ko-','ko-','ko-','ko-','ko-','ko-','ko-'};
a=1.5;
b=1.5;
linewidths={a,a,a,a,a,a,b,a,a};

figurebasename=['figures/SDCRes_N',num2str(N),'_SRange',num2str(100*SRange(1)),'_',num2str(100*SStep),'_',num2str(100*SRange(end))];

%% Make a nice GUI
scrsz = get(0,'ScreenSize');

%%
figure(2)
set(gcf,'Position',[1, 1+scrsz(4)/3, scrsz(3)/2, scrsz(4)/3])
clf
title('S-dependency')


% Plot distribution versus Summed Electrotonic Pathlengths
subplot(1,3,1)
hold on
curves=zeros(length(SRange),Nbin);
for S=1:length(SRange)
    curves(S,:)=accumarray({bin'},P(:,S));
end

for S=1:length(SRange)
    plot(centers,curves(S,:)+(S-1)*0.08,markers{S},'LineWidth',linewidths{S},'MarkerFaceColor','w')
    if(mod(S,2)==1)
        text(centers(1),curves(S,1)+(S-1)*0.08+0.045,['S = ', num2str(SRange(S))],'Fontsize',12)
    end
end
text(centers(1),0.88,['N = ', num2str(N)],'Fontsize',12);
xlabel('Summed Electrotonic Pathlengths','Fontsize',12)
ylabel('Probability Density','Fontsize',12)
set(gca,'YTick',0.08*0:13)
set(gca,'LineWidth',2,'Fontsize',10);
axis auto
axis([centers(1)-0.2,centers(end)+0.1,-0.1,1.1])


% Plot distribution versus Log_2 Multiplicity
subplot(1,3,2)
hold on
offset=0.08;

curves=zeros(length(SRange),N-1);
Mcenters=0:N-2;

for S=1:length(SRange)
    curves(S,:)=accumarray({round(log(10)/log(2)*(M))+1},P(:,S));
end

log2mean=zeros(size(S));

for S=1:length(SRange)
    
    plot(Mcenters,curves(S,:)+(S-1)*offset,markers{S},'LineWidth',linewidths{S},'MarkerFaceColor','w')
    if(mod(S,2)==1)
        text(Mcenters(1),curves(S,1)+(S-1)*offset+0.03,['S = ', num2str(SRange(S))],'Fontsize',12);
    end
end
text(Mcenters(1),0.88,['N = ', num2str(N)],'Fontsize',12);
xlabel('Log_2 Multiplicity','Fontsize',12)
ylabel('Probability Density','Fontsize',12)
set(gca,'YTick',offset*0:13)
%  set(gca,'LineWidth',2);axis auto
set(gca,'LineWidth',2,'Fontsize',10)
axis([Mcenters(1)-1,Mcenters(end)+0.1,-0.1,1])


% Plot distribution versus tree asymmetry index
subplot(1,3,3)
hold on


curves=zeros(length(SRange),NAbin);
for S=1:length(SRange)
    curves(S,:)=accumarray({Abin'},P(:,S));
end

for S=1:length(SRange)
    plot(Acenters,curves(S,:)+(S-1)*0.08,markers{S},'LineWidth',linewidths{S},'MarkerFaceColor','w')
    if(mod(S,2)==1)
        text(Acenters(1),curves(S,1)+(S-1)*0.08+0.03,['S = ', num2str(SRange(S))],'Fontsize',12)
    end
end
text(Acenters(1),0.88,['N = ', num2str(N)],'Fontsize',12);
xlabel('Tree Asymmetry Index','Fontsize',12)
ylabel('Probability Density','Fontsize',12)
set(gca,'YTick',0.08*0:13)
set(gca,'LineWidth',2,'Fontsize',10);
axis auto
axis([0,1,-0.1,1])

%%
figure(3)
set(gcf,'Position',[1, 1+1+scrsz(4)*2/3, scrsz(3)/2, scrsz(4)/3])
clf
hold on

linestyle='.'       % Best for large datasets
linestyle='d-'      % Best for small datasets

subplot(1,3,1)
plot(E(:),linestyle)
ylabel('Summed Electrotonic Pathlengths','Fontsize',12)
xlabel('Tree Index','Fontsize',12)
set(gca,'LineWidth',2,'Fontsize',10)
axis([0.25,double(C)+0.75,min(min(E))*1.05-max(max(E))*0.05,max(max(E))*1.05-min(min(E))*0.05])

subplot(1,3,2)
plot(log(10)/log(2).*M(:),linestyle)
ylabel('Log_2 Multiplicity','Fontsize',12)
xlabel('Tree Index','Fontsize',12)
set(gca,'LineWidth',2,'Fontsize',10)
axis([0.25,C+0.75,0,log(10)/log(2).*max(max(M))+0.5])

subplot(1,3,3)
plot(A(:),linestyle)
ylabel('Tree Asymmetry Index','Fontsize',12)
xlabel('Tree Index','Fontsize',12)
set(gca,'LineWidth',2,'Fontsize',10)
axis([0.25,C+0.75,0,1])
%%
figure(4)
set(gcf,'Position',[1, 1, scrsz(3)/2, scrsz(4)/3])
clf
hold on
subplot(1,3,1)
plot(E(:),A(:),'.')
ylabel('Tree Asymmetry Index','Fontsize',12)
xlabel('Summed Electrotonic Pathlengths','Fontsize',12)
set(gca,'LineWidth',2,'Fontsize',10);axis auto

subplot(1,3,2)
plot(log(10)/log(2).*M(:),E(:),'.')
ylabel('Summed Electrotonic Pathlengths','Fontsize',12)
xlabel('Log_2 Multiplicity','Fontsize',12)
set(gca,'LineWidth',2,'Fontsize',10);
axis auto

subplot(1,3,3)
plot(A(:),log(10)/log(2).*M(:),'.')
ylabel('Log_2 Multiplicity','Fontsize',12)
xlabel('Tree Asymmetry Index','Fontsize',12)
set(gca,'LineWidth',2,'Fontsize',10)
axis auto



%%


format={'png'};
[dummy,size_format]=size(format);
for x=1:1:size_format
    figure(2)
    set(gcf,'PaperPositionMode','auto','PaperOrientation','portrait')
    print(gcf, ['-d',format{x}],'-r120',[figurebasename,'SDependencies_N',num2str(N)])
    figure(3)
    set(gcf,'PaperPositionMode','auto','PaperOrientation','portrait')
    print(gcf, ['-d',format{x}],'-r120',[figurebasename,'TreeAsymmetryMeasureComparison_N',num2str(N)])
    figure(4)
    set(gcf,'PaperPositionMode','auto','PaperSize', [scrsz(3)/2, scrsz(4)/3],'PaperOrientation','portrait')
    print(gcf, ['-d',format{x}],'-r120',[figurebasename,'TreeAsymmetryMeasuresVsIndex_N',num2str(N)])
end

%%
format={'pdf'};
[dummy,size_format]=size(format);
for x=1:1:size_format
    figure(4)
    units= get(gcf,'PaperUnits');
    set(gcf,'PaperUnits','points','PaperSize', [scrsz(3)/2, scrsz(4)/3])
    set(gcf,'PaperPositionMode','auto','PaperSize', [scrsz(3)/2, scrsz(4)/3],'PaperOrientation','portrait')
    print(gcf, ['-d',format{x}],'-r120',[figurebasename,'TreeAsymmetryMeasuresVsIndex_N',num2str(N)])
    set(gcf,'PaperUnits',units)
end



%%
% format={'pdf','jpg'};
% [dummy,size_format]=size(format);
% for x=1:1:size_format
%     figure(2)
%     set(gcf,'PaperPositionMode','auto','PaperOrientation','portrait')
%     saveas(gcf, [figurebasename,'SDependencies_N',num2str(N)],format{x})
%     figure(3)
%     set(gcf,'PaperPositionMode','auto','PaperOrientation','portrait')
%     saveas(gcf, [figurebasename,'TreeAsymmetryMeasureComparison_N',num2str(N)],format{x})
%     figure(4)
%     units= get(gcf,'PaperUnits');
%     set(gcf,'PaperUnits','points')
%     set(gcf,'PaperSize', [scrsz(3)/2, scrsz(4)/3],'PaperPositionMode','auto','PaperOrientation','portrait')
%     saveas(gcf, [figurebasename,'TreeAsymmetryMeasuresVsIndex_N',num2str(N)],format{x})
%     set(gcf,'PaperUnits',units)
% end



