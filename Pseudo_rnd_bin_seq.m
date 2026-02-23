function u=prbs_gen(N, m, a, b)

switch m
    case 3
        coef=[1 3];
     case 4
        coef=[1 4];
     case 5
        coef=[2 5];
     case 6
        coef=[1 6];
     case 7
        coef=[1 7];
     case 8
        coef=[1 2 7 8];
     case 9
        coef=[4 9];
     case 10
        coef=[3 10];
    otherwise
        error('m has to be between 3 and 10');
end


reg=ones(1,m);
u_bits=zeros(1,N);

for k=1:N
   
    u_bits(k)=reg(end);

    fb=0;
    for c=coef
        fb=xor(fb,reg(c));
    end

    reg(2:end)=reg(1:end-1);
    reg(1)=fb;

end

u=a+(b-a)*u_bits;
end


%%
%Generating the signal
Ts=0.01;
N=200;
a=-0.7;
b=0.7;

u_m3=prbs_gen(N,3,a,b);
u_m10=prbs_gen(N,10,a,b);

%validation:step 0.4, aprox 70 samples
u_val=0.4*ones(1,70);

%zeroes between segments
u0_before=zeros(1,150);
u0_between=zeros(1,200);
u0_after=zeros(1,150);

sig=[u0_before,u_m3,u0_between,u_m10,u0_between,u_val,u0_after];

figure(1)
plot(sig);


%Applying signal to DC motor

u=sig(:);

[y_data,~,t_data]=DCMRun.run(u,'Ts',Ts);

save('lab6_motordata.mat','t_data','y_data','u','Ts','sig');

sig
%%
%IDENTIFYING ARX MODELS
load("lab6_motordata.mat")
%index based on signal construction
y_data=detrend(y_data)
idx1_start=150+1;
idx1_end=idx1_start+N-1;

idx2_start=idx1_end+200+1;
idx2_end=idx2_start+N-1;

val_start=idx2_end+200+1;
val_end=val_start+70-1;

%extract id and val data

u_id1=sig(idx1_start:idx1_end).';
u_id2=sig(idx2_start:idx2_end).';
u_val_used=sig(val_start:val_end).';

y_id1=y_data(idx1_start:idx1_end).';
y_id2=y_data(idx2_start:idx2_end).';
y_val_used=y_data(val_start:val_end).';

figure(2)
plot(t_data,y_data);



%build iddata objects


id1=iddata(y_id1,u_id1,Ts);
id2=iddata(y_id2,u_id2,Ts);

val=iddata(y_val_used,u_val_used,Ts);

id1 = detrend(id1);
id2 = detrend(id2);
val = detrend(val);

%na=2; nb=2; nk=1;

model1=arx(id1,[2 2 1]);
model2=arx(id2,[2 2 1]);

figure(3)
compare(val,model1,model2);

PE_m3=2^3-1
PE_m10=2^10-1





