function runs_split=SingleRunPlot_1000ms(target, sigmae, sigmai)
close all

%% Import Data
filename1=sprintf('SpikingNEWNEWv2_p100_sige%1.0f_sigi%1.0f_Trials10_gamma 16.csv', sigmae*10, round(sigmai*10));
filename2=sprintf('SpikingInhibNEWNEWv2_p100_sige%1.0f_sigi%1.0f_Trials10_gamma 16.csv', sigmae*10, round(sigmai*10));
spikingdata=importdata(filename1);
spikingdata_i=importdata(filename2);
S=1;
p=1.00;
trials=10;
Ie=-15.625;
numcells=800;
numcells_i=200;

endtime=2500;

window=100;
start=100;

%% Process Data
runs=cell(S*S,1);
runs_split=cell(S*S,1);
for k=1:S*S
    runs_split{k}=cell(trials,1);
end

markers=find(diff(spikingdata(:,2)));
markers=[0; markers; length(spikingdata(:,1))];
for i=1:length(markers)-1
    runs{i}=spikingdata((1+markers(i)):markers(i+1), :);
end
for i=1:S*S
    temp2=runs{i}(:,3);
    l=length(temp2);
    c=zeros(1,l);
    for j=1:l-1
        c(j)=temp2(j+1)-temp2(j);
    end
    d=find(c<0);
    markers2=[0, d, l];
    for k=1:trials
        runs_split{i}{k}=runs{i}((markers2(k)+1):markers2(k+1),:);
    end
end

%Inhib
runs_i=cell(S*S,1);
runs_split_i=cell(S*S,1);
for k=1:S*S
    runs_split_i{k}=cell(trials,1);
end

markers_i=find(diff(spikingdata_i(:,2)));
markers_i=[0; markers_i; length(spikingdata_i(:,1))];
for i=1:length(markers)-1
    runs{i}=spikingdata_i((1+markers_i(i)):markers_i(i+1), :);
end
for i=1:S*S
    temp2=runs{i}(:,3);
    l=length(temp2);
    c=zeros(1,l);
    for j=1:l-1
        c(j)=temp2(j+1)-temp2(j);
    end
    d=find(c<0);
    markers2_i=[0, d, l];
    for k=1:trials
        runs_split_i{i}{k}=runs{i}((markers2_i(k)+1):markers2_i(k+1),:);
    end
end


%% Synchrony Measure with Shifting Window
window=100;
start=100;
synchrony=zeros(trials, endtime-window-start);
for i=target
    for j=1:trials
        for z=start:endtime-window
            spiketimes=runs_split{i}{j}(:,4);
            spikecell=runs_split{i}{j}(:,3);
            count=1;
            a=spiketimes>z;
            b=spiketimes<z+window;
            spikenum=sum(a.*b);
            formattedspikes=zeros(spikenum, 2);
            for k=1:length(spiketimes)
                if spiketimes(k)>z && spiketimes(k)<(z+window)
                    formattedspikes(count,1)=spiketimes(k)-(z);
                    formattedspikes(count,2)=spikecell(k)+1;
                    count=count+1;
                end
            end
            synchrony(j,z-start+1)=golomb_measure(numcells,formattedspikes,2);
        end
        j
    end
end
% synchrony_smooth=smooth(synchrony,100);

%% Firing Rate with Shifting Window
firingrate=zeros(trials, endtime-window-start);
for i=target
    for j=1:trials
        for z=start:endtime-window
            spiketimes=runs_split{i}{j}(:,4);
            spikecell=runs_split{i}{j}(:,3);
            count=1;
            a=spiketimes>z;
            b=spiketimes<z+window;
            spikenum=sum(a.*b);
            firingrate(j,z-start+1)=(spikenum/numcells)/(window/1000);
        end
        j
    end
end

%% Firing Rate Inhib with Shifting Window
firingrate_i=zeros(trials, endtime-window-start);
for i=target
    for j=1:trials
        for z=start:endtime-window
            spiketimes_i=runs_split_i{i}{j}(:,4);
            spikecell_i=runs_split_i{i}{j}(:,3);
            count=1;
            a=spiketimes_i>z;
            b=spiketimes_i<z+window;
            spikenum_i=sum(a.*b);
            firingrate_i(j,z-start+1)=(spikenum_i/numcells_i)/(window/1000);
        end
        j
    end
end

% %% Load Average Firing Rates
% filename3=sprintf('FiringRateAverage_Window100_p1_Ie-1563_sigmae%1.0f_sigmai%1.0f_ContinuousIncrease.csv', sigmae*1000, (sigmai*1000));
% filename4=sprintf('FiringRateSTD_Window100_p1_Ie-1563_sigmae%1.0f_sigmai%1.0f_ContinuousIncrease.csv', sigmae*1000, (sigmai*1000));
% filename5=sprintf('FiringRateIAverage_Window100_p1_Ie-1563_sigmae%1.0f_sigmai%1.0f_ContinuousIncrease.csv', sigmae*1000, (sigmai*1000));
% filename6=sprintf('FiringRateISTD_Window100_p1_Ie-1563_sigmae%1.0f_sigmai%1.0f_ContinuousIncrease.csv', sigmae*1000, (sigmai*1000));
% 
% filename7=sprintf('SyncAverage_Window100_p1_Ie-1563_sigmae%1.0f_sigmai%1.0f_ContinuousIncrease.csv', sigmae*1000, (sigmai*1000));
% filename8=sprintf('SyncSTD_Window100_p1_Ie-1563_sigmae%1.0f_sigmai%1.0f_ContinuousIncrease.csv', sigmae*1000, (sigmai*1000));
% 
% firingrate_mean=importdata(filename3);
% firingrate_std=importdata(filename4);
% firingrate_mean_i=importdata(filename5);
% firingrate_std_i=importdata(filename6);
% 
% synchrony_mean=importdata(filename7);
% synchrony_std=importdata(filename8);


%% Plot
t=0:.01:endtime;
I=zeros(1, length(t));
for i=1:length(I)
    I(i)=-15.625+((t(i))/(endtime))*2*15.625;
end

grey=[.5 .5 .5];
red=[0.8 0 0];
lightred=[1 .5 .5];
green = [0 0.5 0];

for i=target
    for k=1:trials
        utarget=(1+(i-1)*trials*endtime+(k-1)*endtime):((i-1)*trials*endtime+(k-1)*endtime+endtime);
        
        figure('units','normalized','position',[0 0 .5625 1])

% %         pos1 = [0.15 0.53 0.8 0.07];
        pos1 = [0.15 0.35 0.8 0.6];
        pos2= [0.15 0.2 0.8 0.10];
        
        subplot('Position', pos1)
%         colororder({'r','m'})
        yyaxis left
        scatter(runs_split_i{i}{k}(:,4), runs_split_i{i}{k}(:,3)+numcells, 10, 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', .2);
        hold on
        scatter(runs_split{i}{k}(:,4), runs_split{i}{k}(:,3), 10, 'filled', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', .2);
        hold on
        p1=plot(linspace(0, 2500, 10000), 799.5*ones(10000,1), 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);
        p1.Color(4)=.25;
        hold on
        
        if sigmae==4.4
%             plot(start:1:(endtime-window), (synchrony(k,:)./.45)*1000, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none');
            plot(start:1:(endtime-window), (synchrony(k,:)./.35)*1000, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none');
        elseif sigmae==7.8
%             plot(start:1:(endtime-window), (synchrony(k,:)./.45)*1000, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none');
            plot(start:1:(endtime-window), (synchrony(k,:)./.35)*1000, 'Color', 'b', 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none');
        end
%       axis([0 endtime 0 numcells_i+numcells])
        axis([0 1000 0 numcells_i+numcells])

        set(gca, 'XTick', []);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTick', [0 500 1000]);
%         set(gca, 'YTickLabel', [0 0.225 0.45]);
        set(gca, 'YTickLabel', [0 0.175 0.35]);
        set(gca, 'XTickLabel', []);
        set(gca, 'FontSize', 24);
        ax=gca;
        if sigmae==7.8
            ax.YAxis(1).Color='b';
        elseif sigmae==4.4
            ax.YAxis(1).Color='r';
        end
        hold on
        
        yyaxis right
        plot(start:1:(endtime-window), firingrate(k,:), 'Color', 'k', 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none');
        hold on
        plot(start:1:(endtime-window), firingrate_i(k,:), 'Color', grey, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none');
        hold on        
        set(gca, 'FontSize', 24);
%         axis([0 endtime 0 max([firingrate(k,:), firingrate_i(k,:)])+1])
        axis([0 1000 0 15])
        set(gca, 'XTickLabel', []);
        set(gca, 'FontSize', 24);
        %         legend('Inhibitory Neurons','Excitatory Neurons', 'Location', 'NorthWest')
        %         xlabel('Time, ms', 'FontSize', 34)
        %         ylabel({'Mean Membrane'; 'Voltage Analogue'}, 'FontSize', 30)
        %         plt=gca;
        %         plt.Yaxis(1).Color='r';
        %         plt.Yaxis(2).Color='m';
        ax=gca;
        ax.YAxis(2).Color='k';
     
        subplot('Position', pos2)
        plot(t,I, 'k-', 'LineWidth', 2)
        set(gca, 'FontSize', 24);
        %         xlabel('Time, ms', 'FontSize', 30)
        %         ylabel({'Input to'; 'Excitatory Cells'}, 'FontSize', 24)
%         axis([0 endtime Ie -Ie])      
        axis([0 1000 Ie I(100001)])      
     
        str1=sprintf('SingleRunv2_p%1.0f_Ie%1.0f_sigmae%3.0f_sigmai%3.0f_trial%d_ContinuousIncrease.png', p, Ie*100, runs_split{i}{k}(1,1)*1000, runs_split{i}{k}(1,2)*1000,k);
        %     saveas(gcf, str1)
        set(gcf,'PaperPositionMode','auto')
        print(str1, '-dpng', '-r0');
        
        str2=sprintf('SingleRunv2_p%1.0f_Ie%1.0f_sigmae%3.0f_sigmai%3.0f_trial%d_ContinuousIncrease.eps', p, Ie*100, runs_split{i}{k}(1,1)*1000, runs_split{i}{k}(1,2)*1000,k);
        print(gcf,'-depsc','-painters',str2)

        
        close all
        
    end
end
end