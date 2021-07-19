function HeteroDist(sigmae, sigmai)
close all

%% Import Data
filename1=sprintf('hNEWNEWv2_p100_sige%1.0f_sigi%1.0f_Trials10_gamma 16.csv', sigmae*10, round(sigmai*10));
filename2=sprintf('hInhibNEWNEWv2_p100_sige%1.0f_sigi%1.0f_Trials10_gamma 16.csv', sigmae*10, round(sigmai*10));
h_e=importdata(filename1);
h_i=importdata(filename2);
S=1;
p=1.00;
trials=10;
Ie=-15.625;
numcells=800;
numcells_i=200;

endtime=2500;

window=100;
start=100;

numbin_e=35;
numbin_i=50;

binwidth_e=.75;
binwidth_i=.75;

for i=1:trials
    figure('units','normalized','position',[0 0 .5625 1])
    h_e_now=h_e(1+numcells*(i-1):numcells*i);
    h_i_now=h_i(1+numcells_i*(i-1):numcells_i*i);
    histogram(h_i_now, 'BinWidth', binwidth_i, 'FaceColor', 'k', 'Normalization', 'probability');
    hold on
        if sigmae==4.4
        histogram(h_e_now, 'BinWidth', binwidth_e, 'FaceColor', 'r', 'EdgeColor', 'k', 'Normalization', 'pdf', 'FaceAlpha', .5);
        hold on
    elseif sigmae==7.8
        histogram(h_e_now, 'BinWidth', binwidth_e, 'FaceColor', 'b', 'EdgeColor', 'k', 'Normalization', 'pdf', 'FaceAlpha', .5);
        hold on
    end

    x=linspace(-40, 40, 10000);
    y_e=normpdf(x,0,sigmae);
    y_i=normpdf(x,0,sigmai);
    
    plot(x,y_i, 'LineWidth', 3, 'Color', 'k')
    hold on
    
    if sigmae==4.4
        plot(x,y_e, 'LineWidth', 3, 'Color', 'r')
        hold on
    elseif sigmae==7.8
        plot(x,y_e, 'LineWidth', 3, 'Color', 'b')
        hold on
    end
    
    axis([-40 40 0 .12])
            set(gca, 'FontSize', 30);

    
    str1=sprintf('hDist_sigmae%3.0f_sigmai%3.0f_trial%d_ContinuousIncrease.png', sigmae*1000, sigmai*1000, i);
    %     saveas(gcf, str1)
    set(gcf,'PaperPositionMode','auto')
    print(str1, '-dpng', '-r0');
    
    str2=sprintf('hDist_sigmae%3.0f_sigmai%3.0f_trial%d_ContinuousIncrease.eps', sigmae*1000, sigmai*1000, i);
    print(gcf,'-depsc','-painters',str2)
    
    
    
    
    
    
    
end















end