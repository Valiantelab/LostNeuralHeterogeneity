function [real, imag]=PlotEigenvaluesAndFixedPoints_FIXINGALPHA_080322(filename,K,U, sigmae, sigmai)

%% Import Data
FixedPoints=csvread(filename);

%% Process Data
markers=find(diff(FixedPoints(:,3)));
markers=[0; markers; length(FixedPoints(:,1))];
for i=1:length(markers)-1
    runs{i}=FixedPoints((1+markers(i)):markers(i+1), :);
end

fixed=zeros(length(runs),6);
for i=1:length(runs)
    %% Plot Curves
    close all
    efixed=zeros(0,2);
    ifixed=zeros(0,2);
    
    
    for j=1:length(runs{i}(:,1))
        if runs{i}(j,4)==0
            efixed=[efixed; runs{i}(j,5), runs{i}(j,6)];
        end
        if runs{i}(j,4)==1
            ifixed=[ifixed; runs{i}(j,5), runs{i}(j,6)];
        end
    end
    
    %Find "intersection"
    tol=(U*2)/K;
    count=0;
    check=0;
    for j=1:length(efixed)
        for k=1:length(ifixed)
            dif1=abs(efixed(j,1)-ifixed(k,1));
            dif2=abs(efixed(j,2)-ifixed(k,2));
            %             if dif1<tol && dif2<tol && count==0
            if dif1==0 && dif2==0 && count==0
                fixed(i,1)=efixed(j,1);
                fixed(i,2)=efixed(j,2);
                count=count+1;
            end
            %             if dif1<tol && dif2<tol && count>0
            if dif1==0 && dif2==0 && count>0
                check=0;
                for l=1:count
                    %                     if sqrt((fixed(i,1+(l-1)*2)-efixed(j,1))^2+(fixed(i,2+(l-1)*2)-efixed(j,2))^2)<25*tol
                    %                         check=check+1;
                    %                     end
                    if abs(fixed(i,1+(l-1)*2)-efixed(j,1))<2*tol
                        check=check+1;
                    end
                    if abs(fixed(i,2+(l-1)*2)-efixed(j,2))<2*tol
                        check=check+1;
                    end
                    
                    
                end
                if check==0
                    disp('MULTIPLE FIXED POINTS!!!')
                    disp([i count+1])
                    fixed(i,1+2*count)=efixed(j,1);
                    fixed(i,2+2*count)=efixed(j,2);
                    count=count+1;
                end
            end
        end
    end
end

csvname=sprintf('FixedPointsBifurc_FIXINGALPHA_K%d_U%d_SigE%d_SigI%d.csv', K,U, sigmae*10000, sigmai*10000);
csvwrite(csvname, fixed);

FixedPoints=fixed;









%% Parameters/Storage
a_e=1;
a_i=2;

I=15.625;
U=62;

Gamma=round(1/U,3);
beta=300;

wee=1.6;
wie=-4.7;
wei=3.0;
wii=-0.13;


Isteps=250;



real=zeros(Isteps,8);
imag=zeros(Isteps,8);

%% Loop through fixed points
for i=1:Isteps
    fixedpointnum=ceil(sum(FixedPoints(i,:)~=0)./2);
    
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
        B=-1*(-(a_e+a_i)+(wee*a_e/Gamma)*Fe+(wii*a_i/Gamma)*Fi);
        C=(-a_e+(wee*a_e/Gamma)*Fe)*(-a_i+(wii*a_i/Gamma)*Fi)-((wie*a_i/Gamma)*Fi)*((wei*a_e/Gamma)*Fe);
        
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
%                         if real(i,j)<0 && real(i,j+1)<0
%                             plot(currents(i), real(i,j), 'Color', 'k', 'Marker', '.', 'MarkerSize', 15)
%                             hold on
%                             plot(currents(i), real(i,j+1), 'Color', 'k', 'Marker', '.', 'MarkerSize', 15)
%                             hold on
%                         end
            
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
    
    
    
end
hold on
plot(currents, zeros(1,length(currents)), 'k-', 'LineWidth', 2)
set(gca, 'FontSize', 20);
xlabel('Drive', 'FontSize', 24)
ylabel({'Oscillatory'; 'Frequency (au)'}, 'FontSize', 24)
axis([-I I -22 22])
set(gca, 'XTick', [-15.625, -10.625, -5.625, -0.625, 4.625, 9.625, 14.625])
set(gca, 'XTickLabel', [0 5 10 15 20 25 30]);




%% NOW PLOT FIXED POINT LOCATION TO SHOW BIFURCATION BETTER


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
%                         if imag(i,j)==0 && real(i,j)<0 && real(i,j+1)<0
%                             plot(currents(i), FixedPoints(i,j),  'Color', 'k', 'Marker', '.', 'MarkerSize', 25)
%                             hold on
%                         end
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




str1=sprintf('EigenvaluesBifurcREVISIONS_FIXINGALPHA_SigE%d_SigI%d.png',sigmae*10000, sigmai*10000);
%     saveas(gcf, str1)
set(gcf,'PaperPositionMode','auto')
print(str1, '-dpng', '-r0');

str2=sprintf('EigenvaluesBifurcREVISIONS_FIXINGALPHA_SigE%d_SigI%d.eps',sigmae*10000, sigmai*10000);
print(gcf,'-depsc','-painters',str2)


csvname=sprintf('EigenvaluesBifurcREVISIONSReal_FIXINGALPHA_SigE%d_SigI%d.csv',sigmae*10000, sigmai*10000);
csvwrite(csvname, real);

csvname=sprintf('EigenvaluesBifurcREVISIONSImag_FIXINGALPHA_SigE%d_SigI%d.csv',sigmae*10000, sigmai*10000);
csvwrite(csvname, imag);
end