% plotBEModeCurves

b=1; % basal branching rate
% The parameter E giving the dependence of branch rate on the number of terminal segments will varied
scrsz = get(0,'ScreenSize');
figure(1)
set(gcf,'Position',[1, 1 , scrsz(4) , scrsz(4) ])
clf
hold on

nplotrange=1:5:101;

for E=0:0.25:1
    % Calculate numerical solutions
    [T,Y,mun,sigman,validPoints]=cBEModel(b,E);
    
    % Plot top row: temporal development for the probabilities $p(n,t)$ for
    % $n=1,6,11,...,101$ and three different values of $E$: $E=0$ (left),
    % $E=\frac{1}{2}$ (middle) and $E=1$ (right) in all cases lower
    % $n$-values lead to earlier $p(n,t)$-peaks, while late peaking traces
    % might peak outside the window shown. 
    
    subplot(3,3,1)
    if(E==0)
        plot(T(validPoints),Y(validPoints,nplotrange))
    end
    axis([0,4,0,0.2])
    
    subplot(3,3,2)
    if(E==0.5)
        plot(T(validPoints),Y(validPoints,nplotrange))
    end
    axis([0,25,0,0.2])
    
    subplot(3,3,3)
    if(E==1)
        plot(T(validPoints),Y(validPoints,nplotrange))
    end
    axis([0,25,0,0.2])
    
    
    % Plot middle row: expectation value and variance for the number of terminal
    % segments for three different values of $E$ (matching those in the top
    % row) calculated using the first 1000 $p(n,t)$'s while keeping
    % $p_{high}< 10^{-6}$.  
    subplot(3,3,4)
    hold on
    if(E==0)
        plot(T(validPoints),mun(validPoints),'k-')
        plot(T(validPoints),sigman(validPoints),'k--')
    end
    axis([0,4,0,30])
    
    subplot(3,3,5)
    hold on
    if(E==0.5)
        plot(T(validPoints),mun(validPoints),'k-')
        plot(T(validPoints),sigman(validPoints),'k--')
    end
    axis([0,10,0,30])
    
    subplot(3,3,6)
    hold on
    if(E==1)
        plot(T(validPoints),mun(validPoints),'k-')
        plot(T(validPoints),sigman(validPoints),'k--')
        legend('\mu(n)','\sigma(n )','p_{high}10^7 ','Location','NorthWest')
    end
    axis([0,25,0,30])
    
   
    
    % Plot bottom row left: comparison expectation value (solid grey lines) with
    % mean field prediction (dashed black lines) for five different $E$
    % values: $E= 0,0.25, 0.5, 0.75, 1$; at $E$-values $0, 1$ the mean
    % field solutions coincide with the exact solution.
    subplot(3,2,5)
    hold on
    plot(T(validPoints),mun(validPoints),'k-')
    if(E==0)
        plot(T(validPoints),exp(T(validPoints)),'r--')
    else
        plot(T(validPoints),(E*T(validPoints)+1).^(1/E),'r--')
        legend('\mu(n)','\mu_{MF}(n)','Location','SouthEast')
    end
    axis([0,25,0,30])
    
    
    
    % Plot bottom row right: comparison of mean field solution and numerical
    % results; after an initial growth of the relative error for
    % intermediate values of $E$ , i.e.  $E=  0.25, 0.5, 0.75$, the
    % relative error attenuates. The standard deviation and mean for  $E=0$
    % and $E=1$ correspond within the numerical error with the exact
    % solutions.
    markers=['.';'+';'d';'o';'-'];
   
    subplot(3,2,6)
    hold on
    
    if(E==0 )
        plot(T(validPoints),(exp(T(validPoints))./mun(validPoints)'),markers(1+round(E/0.25)))
    else
        plot(T(validPoints),((E*T(validPoints)+1).^(1/E)./mun(validPoints)'),markers(1+round(E/0.25)))
    end
    
    legend('E=0','E=0.25','E=0.50','E=0.75','E=1','Location','SouthEast')
    axis([0,25,0.99,1.1])
    
end




%% Save figure to file (pdf)

figure(1)
saveas(gcf,'figures/BEModel.pdf')