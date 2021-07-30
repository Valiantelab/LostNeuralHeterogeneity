function [real, imag]=PlotEigenvaluesAndFixedPoints(filename, sigmae, sigmai)


%% Import Data
FixedPoints=csvread(filename);

%% Parameters/Storage
a_e=1;
a_i=2;

I=15.625;
U=62;

% Gamma=0.2/(-I/100);
%or
Gamma=round(1/U,3);
beta=300;

wee=1.6;
wie=-4.7;
wei=3.0;
wii=-0.13;

% sigma_i=1:.5:10.5;
% sigma_e=sigma_i;
% sigma=zeros(S*S,2);

Isteps=250;


% for i=1:length(sigma_e)^2
%     for j=length(FixedPoints(1,:))
%         if FixedPoints(i,j)>10
%             if mod(j,2)==0
%                 FixedPoints(i,j)=0;
%                 FixedPoints(i,j-1)=0;
%             elseif mod(j,2)==1
%                 FixedPoints(i,j)=0;
%                 FixedPoints(i,j+1)=0;
%             end
%         end
%     end
% end

% count=1;
% for i=1:S
%     for j=1:S
%         sigma(count,1)=sigma_e(i);
%         sigma(count,2)=sigma_i(j);
%         count=count+1;
%     end
% end

real=zeros(Isteps,8);
imag=zeros(Isteps,8);

%% Loop through fixed points
for i=1:Isteps
    fixedpointnum=ceil(sum(FixedPoints(i,:)~=0)./2);
%     sigmae=sigma(i,1);
%     sigmai=sigma(i,2);
    
    for j=1:fixedpointnum
        U_e=FixedPoints(i,1+(j-1)*2);
        U_i=FixedPoints(i,j*2);
        
        %% Calculate Convolutions
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
            sum1=sum1+dv*((beta*Gamma*exp(-beta*Gamma*(U_e+v)))/(1+exp(-beta*Gamma*(U_e+v)))^2)...
                *1/sqrt(2*pi*sigmae*sigmae)...
                *exp(-v*v/(2*sigmae*sigmae));
        end
        for k=0:V-1
            v=v_min+k*dv;
            sum2=sum2+dv*((beta*Gamma*exp(-beta*Gamma*(U_i+v)))/(1+exp(-beta*Gamma*(U_i+v)))^2)...
                *1/sqrt(2*pi*sigmai*sigmai)...
                *exp(-v*v/(2*sigmai*sigmai));
        end
        Fe=sum1;
        Fi=sum2;
        
        %% Calculate Eigenvalues
        B=-1*(-2+(wee*a_e/Gamma)*Fe+(wii*a_i/Gamma)*Fi);
        C=(-1+(wee*a_e/Gamma)*Fe)*(-1+(wii*a_i/Gamma)*Fi)-((wie*a_i/Gamma)*Fi)*((wei*a_e/Gamma)*Fe);
        
        cg=(-B+sqrt(B^2-4*C))/2;
        cg0=(-B-sqrt(B^2-4*C))/2;
        cg1=-B/2;
        cg2=sqrt(B^2-4*C)/2;
        cg3=-cg2;
        
        %% Store Eigenvalues
        
        if isreal(cg)==1
            real(i,1+(j-1)*2)=cg;
            real(i,j*2)=cg0;
        else
            real(i,1+(j-1)*2)=cg1;
            real(i,j*2)=cg1;
            imag(i,1+(j-1)*2)=abs(cg2);
            imag(i,j*2)=-abs(cg3);
        end
        
        
    end
    
end

% %% Organize Eigenvalues
% count=1;
% 
% %If there is an imaginary eigenvalue, order it first
% for i=1:Isteps
%     check=find(imag(count,:));
%     if sum(3==check)>0 && sum(1==check)==0
%         real(count, [1 3])=real(count, [3 1]);
%         real(count, [2 4])=real(count, [4 2]);
%         imag(count, [1 3])=imag(count, [3 1]);
%         imag(count, [2 4])=imag(count, [4 2]);
%     end
%     if sum(5==check)>0  && sum(1==check)==0
%         real(count, [1 5])=real(count, [5 1]);
%         real(count, [2 6])=real(count, [6 2]);
%         imag(count, [1 5])=imag(count, [5 1]);
%         imag(count, [2 6])=imag(count, [6 2]);
%     end
% %     if sum(5==check)>0 && sum(1==check)>1 && sum(3==check)==0
% %         real(count, [3 5])=real(count, [5 3]);
% %         real(count, [4 6])=real(count, [6 4]);
% %         imag(count, [3 5])=imag(count, [5 3]);
% %         imag(count, [4 6])=imag(count, [6 4]);
% %     end
%     count=count+1;
% end
% 
% %If there is no imaginary eigenvalue, reorder things so that the first
% %entry is NaN, and also make the second pair the sink
% if length(FixedPoints(1,:))>2
%     count=1;
%     for i=1:Isteps
%         check=find(real(count,:));
%         if imag(count,1)==0 && length(check)==2
%             real(count, [1 3])=real(count, [3 1]);
%             real(count, [2 4])=real(count, [4 2]);
%             imag(count, [1 3])=imag(count, [3 1]);
%             imag(count, [2 4])=imag(count, [4 2]);
%             imag(count, 1)=NaN;
%             imag(count, 2)=NaN;
%             real(count, 1)=NaN;
%             real(count, 2)=NaN;
%         end
%         
%         if imag(count,1)==0 && length(check)==4
%             real(count, [1 3])=real(count, [3 1]);
%             real(count, [2 4])=real(count, [4 2]);
%             imag(count, [1 3])=imag(count, [3 1]);
%             imag(count, [2 4])=imag(count, [4 2]);
%             
%             real(count, [1 5])=real(count, [5 1]);
%             real(count, [2 6])=real(count, [6 2]);
%             imag(count, [1 5])=imag(count, [5 1]);
%             imag(count, [2 6])=imag(count, [6 2]);
%             
%             imag(count, 1)=NaN;
%             imag(count, 2)=NaN;
%             real(count, 1)=NaN;
%             real(count, 2)=NaN;
%         end
% %         count
%         if (real(count,5)<0 && real(count,6)<0) || (real(count,5)>0 && real(count,6)>0)
%             real(count, [4 6])=real(count, [6 4]);
%             real(count, [3 5])=real(count, [5 3]);
%         end
%         
%         if real(count,1)<0 && real(count,2)<0 && imag(count,1)==0
%             real(count, [1 3])=real(count, [3 1]);
%             real(count, [2 4])=real(count, [4 2]);
%             real(count, [1 5])=real(count, [5 1]);
%             real(count, [2 6])=real(count, [6 2]);
%         end
%         
%         count=count+1;
%     end
% end
% 




%% Plot Eigenvalues
currents=linspace(-15.625,15.5,250);

figure('units','normalized','position',[0 0 .5625 1])
green = [0 0.4 0];
blue2= [0.4 0.7 0.9]; 
grey=[.5 .5 .5];
purp=[.4 0 .7];
gold=[.9 .8 0];

subaxis(3,1,2, 'Spacing', 0.01, 'Padding', 0.05, 'Margin', 0.05, 'PaddingBottom', 0.05, 'PaddingTop', 0.00, 'PaddingRight', 0.01, 'PaddingLeft', 0.06);
storereal=zeros(0,0);
count=1;
for i=1:Isteps
    for j=[1 3 5 7]
        if imag(i,j)~=0 && isnan(imag(i,j))==0
            if real(i,j)>0
                plot(currents(i), real(i,j), 'Color', purp, 'Marker', '.', 'MarkerSize', 15)
                hold on
            end
            if real(i,j)<0
                plot(currents(i), real(i,j), 'k.', 'MarkerSize', 15)
                hold on
            end
            storereal(count)=real(i,j);
            count=count+1;
        end
        if real(i,j)~=real(i,j+1)
            
            if real(i,j)<0 && real(i,j+1)<0
                plot(currents(i), real(i,j), 'Color', gold, 'Marker', '.', 'MarkerSize', 15)
                hold on
                plot(currents(i), real(i,j+1), 'Color', gold, 'Marker', '.', 'MarkerSize', 15)
                hold on
            end

            % Change sink to black for panel (d)
%             if real(i,j)<0 && real(i,j+1)<0
%                 plot(currents(i), real(i,j), 'Color', 'k', 'Marker', '.', 'MarkerSize', 15)
%                 hold on
%                 plot(currents(i), real(i,j+1), 'Color', 'k', 'Marker', '.', 'MarkerSize', 15)
%                 hold on
%             end
            
            if real(i,j)*real(i,j+1)<0
                plot(currents(i), real(i,j), 'color', green, 'Marker', '.', 'MarkerSize', 15)
                hold on
                plot(currents(i), real(i,j+1),  'color', green, 'Marker', '.', 'MarkerSize', 10)
                hold on
            end
            
            storereal(count)=real(i,j);
            count=count+1;
            storereal(count)=real(i,j+1);
            count=count+1;
        end
    end
end
hold on
plot(currents, zeros(1,length(currents)), 'k-', 'LineWidth', 2)
set(gca, 'FontSize', 20);
%         xlabel('Time, ms', 'FontSize', 34)
ylabel({'Dampening'; 'Rate'}, 'FontSize', 24)
axis([-I I -3 3])
set(gca, 'XTick', [-15.625, -10.625, -5.625, -0.625, 4.625, 9.625, 14.625])
set(gca, 'XTickLabel', []);








subaxis(3,1,3, 'Spacing', 0.01, 'Padding', 0.05, 'Margin', 0.05, 'PaddingBottom', 0.05, 'PaddingTop', 0.00, 'PaddingRight', 0.01, 'PaddingLeft', 0.06);
for i=1:Isteps
    for j=[1 3 5 7]
        if imag(i,j)~=0 && isnan(imag(i,j))==0
            if real (i,j)>0
                plot(currents(i), imag(i,j), 'Color', purp, 'Marker', '.', 'MarkerSize', 15)
                hold on
                plot(currents(i), imag(i,j+1), 'Color', purp, 'Marker', '.', 'MarkerSize', 15)
                hold on
            end
            if real (i,j)<0
                plot(currents(i), imag(i,j), 'k.', 'MarkerSize', 15)
                hold on
                plot(currents(i), imag(i,j+1), 'k.', 'MarkerSize', 15)
                hold on
            end
        end
    end
%     if imag(i,3)~=0 && isnan(imag(i,1))==0
%         plot(currents(i), imag(i,3), 'b.', 'MarkerSize', 25)
%         hold on
%         plot(currents(i), imag(i,4), 'r.', 'MarkerSize', 25)
%         hold on
%     end
%     if imag(i,5)~=0 && isnan(imag(i,1))==0
%         plot(currents(i), imag(i,5), 'b.', 'MarkerSize', 25)
%         hold on
%         plot(currents(i), imag(i,6), 'r.', 'MarkerSize', 25)
%         hold on
%     end
%     if imag(i,7)~=0 && isnan(imag(i,1))==0
%         plot(currents(i), imag(i,7), 'b.', 'MarkerSize', 25)
%         hold on
%         plot(currents(i), imag(i,8), 'r.', 'MarkerSize', 25)
%         hold on
%     end



end
hold on
plot(currents, zeros(1,length(currents)), 'k-', 'LineWidth', 2)
set(gca, 'FontSize', 20);
xlabel('Drive', 'FontSize', 24)
ylabel({'Oscillatory'; 'Frequency (Hz)'}, 'FontSize', 24)
axis([-I I -22 22])
set(gca, 'XTick', [-15.625, -10.625, -5.625, -0.625, 4.625, 9.625, 14.625])
set(gca, 'XTickLabel', [0 5 10 15 20 25 30]);




%% NOW PLOT FIXED POINT LOCATION TO SHOW BIFURCATION BETTER

% Only for multiple fixed points
% marker=find(FixedPoints(:,3)==0, 1);
% FixedPoints(1:marker-1,[1 5])=FixedPoints(1:marker-1,[5 1]);
% FixedPoints(1:marker-1,[2 6])=FixedPoints(1:marker-1,[6 2]);
% real(1:marker-1,[1 5])=real(1:marker-1,[5 1]);
% real(1:marker-1,[2 6])=real(1:marker-1,[6 2]);
% imag(1:marker-1,[1 5])=imag(1:marker-1,[5 1]);
% imag(1:marker-1,[2 6])=imag(1:marker-1,[6 2]);

subaxis(3,1,1, 'Spacing', 0.01, 'Padding', 0.05, 'Margin', 0.05, 'PaddingBottom', 0.05, 'PaddingTop', 0.00, 'PaddingRight', 0.01, 'PaddingLeft', 0.06);
for i=1:Isteps
    for j=[1 3 5]
        if FixedPoints(i,j)~=0
            if imag(i,j)~=0 && real(i,j)<0
                plot(currents(i), FixedPoints(i,j), 'k.', 'MarkerSize', 25)
                hold on
            end
            if imag(i,j)~=0 && real(i,j)>0
                plot(currents(i), FixedPoints(i,j), 'Color', purp, 'Marker', '.', 'MarkerSize', 25)
                hold on
            end
            if imag(i,j)==0 && real(i,j)*real(i,j+1)<0
                plot(currents(i), FixedPoints(i,j), 'color', green, 'Marker', '.', 'MarkerSize', 25)
                hold on
            end
            
            if imag(i,j)==0 && real(i,j)<0 && real(i,j+1)<0
                plot(currents(i), FixedPoints(i,j),  'Color', gold, 'Marker', '.', 'MarkerSize', 25)
                hold on
            end
            
            % Change sink to black for panel (d)
%             if imag(i,j)==0 && real(i,j)<0 && real(i,j+1)<0
%                 plot(currents(i), FixedPoints(i,j),  'Color', 'k', 'Marker', '.', 'MarkerSize', 25)
%                 hold on
%             end
        end
    end
end
hold on
set(gca, 'FontSize', 20);
% xlabel('Input Current', 'FontSize', 24)
ylabel({'Excitatory'; 'Fixed Point, U_e'}, 'FontSize', 24)
axis([-I I -25 0])
set(gca, 'XTick', [-15.625, -10.625, -5.625, -0.625, 4.625, 9.625, 14.625])
set(gca, 'XTickLabel', []);




str1=sprintf('EigenvaluesBifurcREVISIONS_SigE%d_SigI%d.png',sigmae*10000, sigmai*10000);
%     saveas(gcf, str1)
set(gcf,'PaperPositionMode','auto')
print(str1, '-dpng', '-r0');

str2=sprintf('EigenvaluesBifurcREVISIONS_SigE%d_SigI%d.eps',sigmae*10000, sigmai*10000);
print(gcf,'-depsc','-painters',str2)


csvname=sprintf('EigenvaluesBifurcREVISIONSReal_SigE%d_SigI%d.csv',sigmae*10000, sigmai*10000);
csvwrite(csvname, real);

csvname=sprintf('EigenvaluesBifurcREVISIONSImag_SigE%d_SigI%d.csv',sigmae*10000, sigmai*10000);
csvwrite(csvname, imag);


end
