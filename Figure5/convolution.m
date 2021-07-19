function y=convolution(x,sigmae)
y=zeros(size(x));
%Input scaled x and transform to unscaled Ue for calculation
Ue=(x-472.2227)/35.556;
Gamma=0.016;
gamma=Gamma;
beta=300;

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
    y(i)=sum1;
end

%% Scale Output
y=y*115.2649-.4131;


end