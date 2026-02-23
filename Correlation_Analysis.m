%%obtain data from motor:
Ts=0.01;

% inpu
% t 0
N_zero_init=50;
u_zero_init=zeros(N_zero_init,1);

%random input signal
N_rand=500;
u_rand=(2*0.7)*rand(N_rand,1)-0.7;

%middle zeroes
N_zero_middle=50;
u_zero_middle=zeros(N_zero_middle,1);

%final signal
N_step=70;
step_magnitude=0.2;
u_step=step_magnitude*ones(N_step,1);

u=[u_zero_init;u_rand;u_zero_middle;u_step];
N_total=length(u);
t_final=(N_total-1)*Ts;
t=(0:Ts:t_final)';

[y_data,~,t_data]=DCMRun.run(u,'Ts',Ts);
save('lab4_motordata.mat','t_data','y_data','u','Ts');




%%
%main:
load("lab4_motordata.mat")
figure(1)
plot(t_data,u);
figure(2)
plot(t_data,y_data);

N_zero_init=50;
N_rand=500;
N_zero_middle=50;
N_step=70;
N_total=670;
N_id=N_rand;%id set length

%identification
start_id=N_zero_init+1;
end_id=N_zero_init+N_rand;
u_id=u(start_id:end_id);
y_id=y_data(start_id:end_id);


%validation
start_val=N_zero_init+N_rand+N_zero_middle+1;
end_val=N_total;
u_val=u(start_val:end_val);
y_val=y_data(start_val:end_val);

%identification zero_mean or not
figure(3)
plot(u_id);hold on
plot(y_id,'r');

%time vector
t_id=t_data(start_id:end_id);

%check if data is zero-mean or not
mean_u_id=mean(u_id);
mean_y_id=mean(y_id);
mean_u_id
mean_y_id

%get zero-mean signals
u_id_zm=detrend(u_id,0);
y_id_zm=detrend(y_id,0);

figure(4);
plot(t_id,u_id_zm,'g');hold on;
plot(t_id,y_id_zm,'b');

%now check for average to be around 0
mean_u_id_check=mean(u_id_zm)
mean_y_id_check=mean(y_id_zm)

N_id=length(u_id_zm);
M_max=50;

r_u=zeros(M_max+1,1);
r_yu=zeros(M_max+1,1);



u_id_zm=u_id_zm(:);
y_id_zm=y_id_zm(:);


for tau=0:M_max
    length_sum=N_id-tau;
    Normalizare=1/length_sum;

    u_shifted=u_id_zm(tau+1:N_id); %u(k+tau)
    u_base=u_id_zm(1:length_sum);

    y_shifted=y_id_zm(tau+1:N_id);

    %r_u(tau)
    sum_u=sum(u_shifted.*u_base);
    r_u(tau+1)=sum_u*Normalizare;

    %r_yu(tau)
    sum_yu=sum(y_shifted.*u_base);
    r_yu(tau+1)=sum_yu*Normalizare;

end
figure(5)
plot(0:M_max,r_u,'r')
figure(6)
plot(0:M_max,r_yu);

%solving syst of lin eq to obtain FIR model
M=21;
T=min(M+15,length(r_u));
if T<=M
    T=min(M+1, length(r_u));
end

ru=r_u(:);
ryu=r_yu(:);

%Ru matrix
Ru=zeros(T,M);
for i=1:T
    for j=1:M
        Ru(i,j)=ru(abs(i-j)+1);
    end
end


%free erms vector:
Ryu=ryu(1:T);

h=Ru\Ryu;

figure(11)
plot(0:M-1,h);

%convolution to simulate model response
Ts=0.01;
t_eq=(0:length(u)-1).'*Ts;

u_all_zm=u-mean_u_id;
y_hat_zm=conv(u_all_zm,h);
y_hat=y_hat_zm(1:length(u))+mean_y_id;
%val and id 
y_hat_id=y_hat(start_id:end_id);
y_hat_val=y_hat(start_val:end_val);

fit_id=1-norm(y_id-y_hat_id)^2/norm(y_id-mean(y_id))^2
fit_val=1-norm(y_val-y_hat_val)^2/norm(y_val-mean(y_val))^2

rmse_id=sqrt(mean((y_id-y_hat_id).^2))
rmse_val=sqrt(mean((y_val-y_hat_val).^2))

figure(7)
plot(t_data,y_data,'b');hold on
plot(t_eq,y_hat,'r');

figure(8)
plot(t_data(start_id:end_id),y_id,'b');hold on
plot(t_eq(start_id:end_id),y_hat_id,'r');

figure(9)
plot(t_data(start_val:end_val),y_val,'b');hold on
plot(t_eq(start_val:end_val),y_hat_val,'r');

M_max=min(50, length(r_u)-1);
M_grid=1:1:M_max;

u_all_zm=u-mean_u_id;

fit_id_grid=zeros(length(M_grid),1);
fit_val_grid=zeros(length(M_grid),1);
rmse_id_grid=zeros(length(M_grid),1);
rmse_val_grid=zeros(length(M_grid),1);

u_all_zm=u_all_zm(:);
y_id=y_id(:);
y_val=y_val(:);

for k=1:length(M_grid)
    M=M_grid(k);
    T=M+10;
    
    if T>length(ru)
        T=length(ru);
    end

    %Ru matrix
Ru_current=zeros(T,M);
for i=1:T
    for j=1:M
        Ru_current(i,j)=ru(abs(i-j)+1);
    end
end

    %free term vector
Ryu_current=ryu(1:T);

h=Ru_current\Ryu_current;

y_hat_zm=conv(u_all_zm,h);
y_hat=y_hat_zm(1:N_total)+mean_y_id;

%separate into id and val
y_hat_id=y_hat(start_id:end_id);
y_hat_val=y_hat(start_val:end_val);

fit_id_grid(k)=1-norm(y_id-y_hat_id)^2/norm(y_id-mean(y_id))^2;
fit_val_grid(k)=1-norm(y_val-y_hat_val)^2/norm(y_val-mean(y_val))^2;

rmse_id_grid(k)=sqrt(mean((y_id-y_hat_id).^2))
rmse_val_grid(k)=sqrt(mean((y_val-y_hat_val).^2))


end
figure(10)
plot(M_grid, fit_id_grid,'r');hold on
plot(M_grid, fit_val_grid,'b');
plot(M_grid, rmse_id_grid,'g');hold on
plot(M_grid, rmse_val_grid,'y');

M_large=min(50,length(ru)-1);
L_rows=10;
T_large=min(M_large+L_rows,length(ru));

RuL=zeros(T_large,M_large);
for i=1:T_large
    for j=1:M_large
        RuL(i,j)=ru(abs(i-j)+1);
    end

end

h_large=RuL\ryu(1:T_large);
hmax=max(abs(h_large));
prag=0.02*hmax;
M_trans=find(abs(h_large)>prag,1,'last');
if isempty(M_trans)
    M_trans=1;
end

bestfit=max(fit_val_grid);
idx=find(fit_val_grid>=0.98*bestfit,1,'first');
M_gen=M_grid(idx);

M_opt=max(M_trans,M_gen);


