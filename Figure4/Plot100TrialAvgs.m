function runs_split=Plot100TrialAvgs(target, filename1, filename2)
close all

%% Import Data
spikingdata=importdata(filename1);
spikingdata_i=importdata(filename2);
S=1;
p=1.00;
trials=100;
Ie=-15.625;
numcells=800;
numcells_i=200;
endtime=2500;


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


%% Mean and STD
len=length(synchrony);
synchrony_mean=zeros(1,len);
synchrony_std=zeros(1,len);
for i=1:len
    synchrony_temp=synchrony(:,i);
    synchrony_temp=synchrony_temp(~isnan(synchrony_temp));
    synchrony_mean(i)=mean(synchrony_temp);
    synchrony_std(i)=std(synchrony_temp);
end

% synchrony_mean=mean(synchrony);
firingrate_mean=mean(firingrate);
% synchrony_std=std(synchrony);
firingrate_std=std(firingrate);
% synchrony_mean(isnan(synchrony_mean))=0;
% synchrony_std(isnan(synchrony_std))=0;
firingrate_mean_i=mean(firingrate_i);
firingrate_std_i=std(firingrate_i);

%% Plot
t=0:.01:endtime;
I=zeros(1, length(t));
for i=1:length(I)
    I(i)=-15.625+((t(i))/(endtime))*2*15.625;
end

figure('units','normalized','position',[0 0 1 1])
pos1 = [0.15 0.45 0.8 0.2];
pos3 = [0.15 0.225 0.8 0.2];
%         pos2 = [0.15 0.225 0.8 0.10];

pos4 = [0.15 0.10 0.8 0.10];

subplot('Position', pos1)
plot(start:1:(endtime-window), synchrony_mean, 'b-', 'LineWidth', 3);
hold on
p1=plot(start:1:(endtime-window), synchrony_mean-synchrony_std, 'b-', 'LineWidth', 2);
hold on
p2=plot(start:1:(endtime-window), synchrony_mean+synchrony_std, 'b-', 'LineWidth', 2);
p1.Color(4)=0.25;
p2.Color(4)=0.25;
set(gca, 'FontSize', 24);
axis([0 endtime 0 0.45])
set(gca, 'XTickLabel', []);

subplot('Position', pos3)
plot(start:1:(endtime-window), firingrate_mean, 'm-', 'LineWidth', 3);
hold on
p1=plot(start:1:(endtime-window), firingrate_mean-firingrate_std, 'm-', 'LineWidth', 2);
hold on
p2=plot(start:1:(endtime-window), firingrate_mean+firingrate_std, 'm-', 'LineWidth', 2);
p1.Color(4)=0.25;
p2.Color(4)=0.25;
set(gca, 'FontSize', 24);
axis([0 endtime 0 max(firingrate_mean+firingrate_std)+1])
set(gca, 'XTickLabel', []);


%         subplot('Position', pos2)
%         plot(start:1:(endtime-window), synchrony(k,:), 'LineWidth', 2);
%         set(gca, 'FontSize', 24)
% %         ylabel({'Synchrony Measure'; '100 ms Window'}, 'FontSize', 24)
%         axis([0 endtime 0 .45])
%         set(gca, 'XTickLabel', []);


subplot('Position', pos4)
plot(t,I, 'k-', 'LineWidth', 2)
set(gca, 'FontSize', 24);
%         xlabel('Time, ms', 'FontSize', 30)
%         ylabel({'Input to'; 'Excitatory Cells'}, 'FontSize', 24)
axis([0 endtime Ie -Ie])

str1=sprintf('AvgMeasures_Window%1.0d_p%1.0f_Ie%1.0f_sigmae%3.0f_sigmai%3.0f_ContinuousIncrease.png', window, p, Ie*100, runs_split{1}{1}(1,1)*1000, runs_split{1}{1}(1,2)*1000);
%     saveas(gcf, str1)
set(gcf,'PaperPositionMode','auto')
print(str1, '-dpng', '-r0');

str2=sprintf('SyncAverage_Window%1.0d_p%1.0f_Ie%1.0f_sigmae%3.0f_sigmai%3.0f_ContinuousIncrease.csv', window, p, Ie*100, runs_split{1}{1}(1,1)*1000, runs_split{1}{1}(1,2)*1000);
csvwrite(str2, synchrony_mean);
str3=sprintf('SyncSTD_Window%1.0d_p%1.0f_Ie%1.0f_sigmae%3.0f_sigmai%3.0f_ContinuousIncrease.csv', window, p, Ie*100, runs_split{1}{1}(1,1)*1000, runs_split{1}{1}(1,2)*1000);
csvwrite(str3, synchrony_std);
str4=sprintf('FiringRateAverage_Window%1.0d_p%1.0f_Ie%1.0f_sigmae%3.0f_sigmai%3.0f_ContinuousIncrease.csv', window, p, Ie*100, runs_split{1}{1}(1,1)*1000, runs_split{1}{1}(1,2)*1000);
csvwrite(str4, firingrate_mean);
str5=sprintf('FiringRateSTD_Window%1.0d_p%1.0f_Ie%1.0f_sigmae%3.0f_sigmai%3.0f_ContinuousIncrease.csv', window, p, Ie*100, runs_split{1}{1}(1,1)*1000, runs_split{1}{1}(1,2)*1000);
csvwrite(str5, firingrate_std);
str6=sprintf('FiringRateIAverage_Window%1.0d_p%1.0f_Ie%1.0f_sigmae%3.0f_sigmai%3.0f_ContinuousIncrease.csv', window, p, Ie*100, runs_split{1}{1}(1,1)*1000, runs_split{1}{1}(1,2)*1000);
csvwrite(str6, firingrate_mean_i);
str7=sprintf('FiringRateISTD_Window%1.0d_p%1.0f_Ie%1.0f_sigmae%3.0f_sigmai%3.0f_ContinuousIncrease.csv', window, p, Ie*100, runs_split{1}{1}(1,1)*1000, runs_split{1}{1}(1,2)*1000);
csvwrite(str7, firingrate_std_i);

close all

end