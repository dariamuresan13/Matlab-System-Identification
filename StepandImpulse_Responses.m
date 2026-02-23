%date_motor:
Ts = 0.01;

t1 = 0.3; 
t2 = 0.7; 
t3 = 0.3;
t4 = 0.7; 
t5 = 0.3; 
t6 = 0.7;
t7 = 0.3; 

t_final = t1 + t2 + t3 + t4 + t5 + t6 + t7;

t = (0:Ts:t_final-Ts);

u = zeros(size(t));

start1 = t1;
stop1  = start1 + t2;
start2 = stop1 + t3;
stop2  = start2 + t4;
start3 = stop2 + t5;
stop3  = start3 + t6;

u(t >= start1 & t < stop1) = 0.7;
u(t >= start2 & t < stop2) = 0.4;
u(t >= start3 & t < stop3) = -0.5;

[y_data, ~, t_data] = DCMRun.run(u,'Ts',Ts);

save('motor_data.mat', 't_data', 'y_data', 'u', 'Ts');

%%
%impuls:
Ts=0.01;
t=0:Ts:5.0;
u_const=0.1*ones(size(t));
[y_data_const, ~, t_data_const] = DCMRun.run(u_const,'Ts',Ts);
figure
plot(t_data_const,y_data_const);


t_ss=0.3; %t for ss regime 
n_ss=round(t_ss/Ts); %samples for ss
n_imp=2; %impulse time(2 samples =0.02s)

u_seg_ss=0.1*ones(n_ss,1); %segment ss
u_seg_imp=1.0*ones(n_imp,1); %segment impulse

u_impuls=[
    u_seg_ss;
     u_seg_imp;
      u_seg_ss;
       u_seg_imp;
        u_seg_ss;
         u_seg_imp;
          u_seg_ss;

];

N=length(u_impuls);
t=(0:N-1)'*Ts;

[y_data_imp, ~, t_data_imp] = DCMRun.run(u_impuls,'Ts',Ts);

save('motor_data_impulse.mat', 't_data_imp', 'y_data_imp', 'u_impuls', 'Ts');



%%
%Main code:

load("motor_data.mat")
load("motor_data_impulse.mat")

figure(1)
plot(t_data,u, 'g');

figure(2)
plot(t_data,y_data, 'r');

%FIRST PART
dt=diff(t_data);

delta_dt=max(dt)-min(dt);

delta_dt

if delta_dt>0.01
    N = length(t_data);
    t_uniform=(0:N-1)' * Ts;
    y_data_resampled=interp1(t_data,y_data,t_uniform);
    u_resampled=interp1(t_data,u,t_uniform);
    t_data=t_uniform;
    y_data=y_data_resampled;
    u=u_resampled;
end


%id data from first step
identification=(t_data>0.3 & t_data<1.0);
t_id=t_data(identification);
u_id=u(identification);
y_id=y_data(identification);


%transfer function model with first step signal and response
u_ss=0.7;
t_start_step=0.53;
u_0=0.0;

y_0=mean(y_data(t_data<t_start_step));
y_ss=mean(y_data(t_data>0.8 & t_data<=1.0));


%gain
K=(y_ss-y_0)/(u_ss-u_0);

%time constant
y_T=y_0+0.632*(y_ss-y_0);

ix_T = find(y_id >= y_T, 1, 'first');

T = t_id(ix_T) -t_start_step;


%TF:
H=tf(K,[T 1]);

K
T
H

%Validating data(steps 2 and 3)


inx_start_val = find(t_data<1, 1, 'last') + 1;
inx_stop_val = length(t_data);
validation = inx_start_val:inx_stop_val;
t_val=t_data(validation);
u_val=u(validation);
y_val=y_data(validation);

%applying model to validation data
Ts=0.01;
Nt = length(t_data);

t_sim_vector = (0:Nt-1)' * Ts;

y_model_full = lsim(H, u, t_sim_vector);
y_model = y_model_full(validation);

figure(3)
plot(t_val,y_val,'b')
hold on
plot(t_val,y_model,'r')

MSE=mean((y_val-y_model).^2);
MSE
%% 

%SECOND PART


Ts=0.01;
t_ss_duration=0.3; 
n_ss=round(t_ss_duration/Ts); %samples for ss
n_imp=2;
u_ss_imp=0.1;

dt_imp=diff(t_data_imp);

delta_dti=max(dt_imp)-min(dt_imp);

delta_dti

if delta_dti>0.01
    Nimp = length(t_data_imp);
    t_uniform_imp=(0:Nimp-1)' * Ts;
    y_data_imp=interp1(t_data_imp,y_data_imp,t_uniform_imp);
    u_impuls=interp1(t_data_imp,u_impuls,t_uniform_imp);
    t_data_imp=t_uniform_imp;
   
end

%ideal dirac impulse area=1
%our impulse has amplitude of 1
%delta t=n_imp*Ts=2*0.01=0.02
%A*delta t=1*0.02=0.02
%for area=1:
%Anew*delta t=1 so Anew=1/delta t=1/0.02=50
%alpha=Aideal/A=50/1.0=50

alpha = 50;

%transfer f from 1 st impulse
y_0_imp=mean(y_data_imp(1:n_ss));%ss

idx_start_id=n_ss+1;
idx_end_id=n_ss+n_imp+n_ss; %first impulse

t_id_imp=t_data_imp(idx_start_id:idx_end_id);
y_id_imp=y_data_imp(idx_start_id:idx_end_id);

y_id_norm=y_id_imp-y_0_imp;



%impulse peak
[y_max_imp, idx_max_imp]=max(y_id_norm);
t_peak_imp=t_id_imp(idx_max_imp);%t1

y_T_target=0.368*y_max_imp; 

%output after peak (going down)
y_decay=y_id_norm(idx_max_imp:end);
t_decay=t_id_imp(idx_max_imp:end);

idx_T_imp=find(y_decay<=y_T_target,1,'first');
t_T_imp_decay=t_decay(idx_T_imp);%t2

T_imp=t_T_imp_decay-t_peak_imp;

%K and TF for impulse
u_ss_val=0.1;
K_imp=y_0_imp/u_ss_val;

H_imp=tf(K_imp,[T_imp 1]);
H_imp

%validation data
idx_start_val_imp=idx_end_id+1;
idx_end_val_imp=length(u_impuls);

u_val_imp=u_impuls(idx_start_val_imp:idx_end_val_imp);
t_val_imp=t_data_imp(idx_start_val_imp:idx_end_val_imp);
y_val_imp=y_data_imp(idx_start_val_imp:idx_end_val_imp);

%state space

A=-1/T_imp;
B=K_imp/T_imp;
C=1;
D=0;
sys_ss=ss(A,B,C,D);
x0=y_data_imp(idx_end_id);%initial state
y_model_imp=lsim(sys_ss,u_val_imp,t_val_imp,x0);
figure(4)
plot(t_val_imp,y_val_imp,'b');hold on
plot(t_val_imp,y_model_imp,'r');

figure
plot(u_impuls);
figure
plot(y_val_imp);