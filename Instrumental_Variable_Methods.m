clear;clc;close all;
load("iv_3.mat");

figure(1)
plot(id);
title('identification data')
figure(2)
plot(val);
title('validation data')


u_id=id.InputData(:);
y_id=id.OutputData(:);

u_val=val.InputData(:);
y_val=val.OutputData(:);

na=n;
nb=n; 
nk=1; 

N=length(y_id);
max_order=max(na,nb);
Phi=zeros(N-max_order,na+nb);
Y=zeros(N-max_order,1);


%building regressor matrix phi

for k=max_order+1:N
    row=[ ];
    for i=1:na
        row=[row,-y_id(k-i)];
    end

    for j=1:nb
        row=[row,u_id(k-j)];
    end

    Phi(k-max_order, :)=row;
    Y(k-max_order)=y_id(k);
end

%parameters matrix theta
     theta_ARX=Phi\Y

%we separate coefficients
     a=theta_ARX(1:na);
     b=theta_ARX(na+1:end);


%ARX Simulation     
y_sim=zeros(N,1);

for k=1:N
    Asum=0;
    for i=1:na
        if k-i>0
            Asum=Asum-a(i)*y_sim(k-i);
        end
    end

    Bsum=0;
    for j=1:nb
        if k-j>0
            Bsum=Bsum+b(j)*u_id(k-j);
        end
    end

    y_sim(k)=Asum+Bsum;
end

Z=zeros(N-max_order, na+nb);

for k=max_order+1:N
    row=[];
    for i=1:na
        row=[row, -y_sim(k-i)]; 
    end
    for j=1:nb
        row=[row, u_id(k-j)]; 
    end

    Z(k-max_order,:)=row;
end

% IV Ecuations

Phi_tilde=(Z'*Phi)/(N-max_order);
Y_tilde=(Z'*Y)/(N-max_order);

theta_IV=Phi_tilde\Y_tilde; 


a_iv=theta_IV(1:na); 
b_iv=theta_IV(na+1:end);


Ts=val.Ts;
A_poly=[1;a_iv]'; 
B_poly=[0;b_iv]';

model_iv=idpoly(A_poly,B_poly,[],[],[],0,Ts);

arx_model=arx(id,[na nb nk]);
figure
compare(val,model_iv,arx_model);

