clear;
clc;
close all;
load("arxdata_4.mat");
figure(1)
plot(id);
title('identification data')
figure(2)
plot(val);
title('validation data')

%get identification and val data
u_id=id.InputData(:);
y_id=id.OutputData(:);

u_val=val.InputData(:);
y_val=val.OutputData(:);



%identify based on step response
na=4;
nb=3;
nk=1; %delay 0



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
     theta=Phi\Y

%we separate autoregressive and exogenous coefficients
     a=theta(1:na);
     b=theta(na+1:end);


     
Nval=length(y_val);
y_hat=zeros(Nval,1);

for k=1:Nval
    Asum=0;
    for i=1:na
        if k-i>0
            Asum=Asum-a(i)*y_hat(k-i);
        end
    end

    Bsum=0;
    for j=1:nb
        if k-j>0
            Bsum=Bsum+b(j)*u_val(k-j);
        end
    end

    y_hat(k)=Asum+Bsum;
end



Ts=val.Ts;
sim=iddata(y_hat,u_val,Ts);
figure(5)
compare(sim,val)

% pred
y_pred=zeros(Nval,1);
for k=1:Nval
    Asump=0;
    for i=1:na
        if k-i>0
            Asump=Asump-a(i)*y_val(k-i);
        end
    end
    Bsump=0;
    for j=1:nb
         if k-j>0
             Bsump=Bsump+b(j)*u_val(k-j);
         end
    end
    y_pred(k)=Asump+Bsump;
end




pre=iddata(y_pred,u_val,Ts);
figure(6)
compare(pre,val)

figure(7)
model=arx(id,[na,nb,nk]);
compare(model,val);
