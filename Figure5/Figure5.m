close all

sigmae=4.4;
sigmae2=7.8;

a_e=1;
a_i=2;

Gamma=0.016;
gamma=Gamma;
beta=300;

wee=1.6;
wie=-4.7;
wei=3.0;
wii=-0.13;

Ue=-15.625:.625:15.625;
Fe=zeros(1,length(Ue));
Fe2=zeros(1,length(Ue));

%% Calculate Convolutions
for i=1:length(Ue)
    U_e=Ue(i);
    v_min0=-1;
    v_max0=1;
    v_min=v_min0/Gamma;
    v_max=v_max0/Gamma;
    V=1000;
    dv=abs(v_max-v_min)/V;
    sum1=0;
    sum2=0;
    for k=0:V-1
        v=v_min+k*dv;
        sum1=sum1+dv*(1/(1+exp(-beta*gamma*(U_e+v))))...
            *1/sqrt(2*pi*sigmae*sigmae)...
            *exp(-v*v/(2*sigmae*sigmae));
    end
    for k=0:V-1
        v=v_min+k*dv;
        sum2=sum2+dv*(1/(1+exp(-beta*gamma*(U_e+v))))...
            *1/sqrt(2*pi*sigmae2*sigmae2)...
            *exp(-v*v/(2*sigmae2*sigmae2));
    end
    Fe(i)=sum1;
    Fe2(i)=sum2;
    
end

%% Extract Example Sigmoids
Ue2=-20:.1:20;
Fe_ex=zeros(100, length(Ue2));
Fe2_ex=zeros(100,length(Ue2));
for i=1:100
    h_1=normrnd(0,1)*sigmae;
    h_2=normrnd(0,1)*sigmae2;
    Fe_ex(i,:)=(1./(1+exp(-beta.*gamma.*(Ue2-h_1))));
    Fe2_ex(i,:)=(1./(1+exp(-beta.*gamma.*(Ue2-h_2))));
end


red=[0.8 0 0];
blue=[0 0.1 0.4];

figure('units','normalized','position',[0 0 1 1])
for i=1:100
    p1=plot(Ue2,Fe_ex(i,:), 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2, 'Marker', 'none');
    hold on
    p2=plot(Ue2,Fe2_ex(i,:), 'Color', 'b', 'LineStyle', '-', 'LineWidth', 2, 'Marker', 'none');
    hold on
    p1.Color(4)=0.3;
    p2.Color(4)=0.3;
end
p3=plot(Ue, Fe, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 4, 'Marker', 'none');
hold on
p4=plot(Ue, Fe2, 'Color', 'b', 'LineStyle', '-', 'LineWidth', 4, 'Marker', 'none');
hold on
set(gca,'box', 'off')
set(gca,'TickDir','out')
 
xlabel('Membrane Potential Analogue', 'FontSize', 30)
ylabel('Firing Probability', 'FontSize', 30)
axis([-20 20 0 1])
set(gca, 'FontSize', 24);
legend([p3 p4 p1 p2], {'Average Epileptic', 'Average Non-epileptic', 'Single Neurons, Epileptic Range', 'Single Neurons, Non-epileptic range'}, 'FontSize', 20, 'Location', 'SouthEast', 'NumColumns', 2)

str1=sprintf('PopulationFICurves_v2.png');
%     saveas(gcf, str1)
set(gcf,'PaperPositionMode','auto')
print(str1, '-dpng', '-r0');

str2=sprintf('PopulationFICurves_v2.eps');
print(gcf,'-depsc','-painters',str2)


%% Plot with Experimental Data
figure('units','normalized','position',[0 0 1 1])
currents=[50 100 150 200 250];
freq_e=[0 .727 5.100 8.136 10.955];
freq_n=[5.410 9.682 13.890 19.414 23.924];
std_e=[0 0.331687396 1.821487903 1.975958811 2.235228833];
std_n=[1.604819909 1.317792496 1.503691376 1.877708134 2.474415305];


x=currents;

blue2=[.1 .5 1];
red2=[1 .4 .2];
green=[0.3 0.4 0];
x1=linspace(currents(1)-1,currents(end)+1,500);
x1_s=x1*0.0278-13.2372;
fit_e=convolution(x1,5.039);
fit_n=convolution(x1,7.771);
fit_e_s=0.008844674803545*fit_e+.003584;
fit_n_s=0.008844674803545*fit_n+.003584;

p1=errorbar(currents, freq_e, std_e, 'Color', 'r', 'Marker', '.', 'MarkerSize', 40, 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 10);
hold on
p2=errorbar(currents, freq_n, std_n, 'Color', 'b', 'Marker', '.', 'MarkerSize', 40, 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 10);
hold on
axis([49 251 0 27])
set(gca,'box', 'off')
set(gca,'TickDir','out')
ylabel('Average Firing Frequency (Hz)', 'FontSize', 30)
xlabel('Input Current (pA)', 'FontSize', 30)
ax1=gca;

ax1_pos=ax1.Position;
ax2=axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
p3=line(Ue, Fe, 'Parent', ax2, 'Color', 'r', 'Marker', 'none', 'LineStyle', '-', 'LineWidth', 3);
hold on
p4=line(Ue, Fe2, 'Parent', ax2, 'Color', 'b', 'Marker', 'none', 'LineStyle', '-', 'LineWidth', 3);
hold on
p5=line(x1_s, fit_e_s, 'Parent', ax2, 'Color', red2, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', 3);
hold on
p6=line(x1_s, fit_n_s, 'Parent', ax2, 'Color', blue2, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', 3);
ax2.XColor=green;
ax2.YColor=green;
xlabel('Membrane Potential Analogue', 'FontSize', 30)
ylabel('Firing Probability', 'FontSize', 30)

axis([-11.875 -6.25 0.0 0.235929511111111])
set(gca,'TickDir','out')
set(ax1, 'FontSize', 24);
set(ax2, 'FontSize', 24);
legend([p1 p2 p3 p4 p5 p6], {'Epileptic, Experiment', 'Non-epileptic, Experiment', ...
    'Epileptic Heterogeneity, Experiment', 'Non-epileptic Heterogeneity, Experiment', ...
    'Epileptic Best Fit', 'Non-epileptic Best Fit'}, 'FontSize', 20, 'Location', 'NorthWest', 'NumColumns', 3)

str1=sprintf('PopulationFICurves_DataAndFit_v2.png');
%     saveas(gcf, str1)
set(gcf,'PaperPositionMode','auto')
print(str1, '-dpng', '-r0');

str2=sprintf('PopulationFICurves_DataAndFit_v2.eps');
print(gcf,'-depsc','-painters',str2)
