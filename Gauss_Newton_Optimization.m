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


%Generating the signal
Ts=0.01;
N=200;
u_prbs=prbs_gen(N,10,-0.7,0.7);
u_step=0.4*ones(1,70);


%zeroes between segments
u0_before=zeros(1,150);
u0_between=zeros(1,200);


sig=[u0_before,u_prbs,u0_between,u_step];

%Applying signal to DC motor

u=sig(:);

[y_data,~,t_data]=DCMRun.run(u,'Ts',Ts);

save('lab7_motordata.mat','t_data','y_data','u','Ts','sig');

%%

load('lab7_motordata.mat')
figure(1)
plot(sig);

figure(2)
plot(t_data,y_data)

%identification data

N=200;
start_prbs=150+1;
end_prbs=150+N;

u_id=u(start_prbs:end_prbs);
y_id=y_data(start_prbs:end_prbs);

su=std(u_id);
sy=std(y_id);

y_id= y_id /sy;
u_id = u_id /su;



%Algorithm
alpha=0.1;
delta=1e-5;
lmax=100;
nk=4;

f0=0.1;
b0=0.1;

theta=[f0,b0];

N=length(y_id);

%gauss newton algorithm


for l=1:lmax
    f=theta(1);
    b=theta(2);

    epsk=zeros(N,1);
    dfk=zeros(N,1);
    dbk=zeros(N,1);


    for k=nk+1:N
        epsk(k)=y_id(k)+f*y_id(k-1)-b*u_id(k-nk)-f*epsk(k-1);

        dfk(k)=y_id(k-1)-epsk(k-1)-f*dfk(k-1);

        dbk(k)=-u_id(k-nk)-f*dbk(k-1);

    end

    idx=nk+1:N;


    %gradient
    grad=(2/length(idx))*[sum(epsk(idx).*dfk(idx)); sum(epsk(idx).*dbk(idx))];
    %approximate Hessian
    H=(2/length(idx))*[sum(dfk(idx).^2),sum(dfk(idx).*dbk(idx)); sum(dfk(idx).*dbk(idx)),sum(dbk(idx).^2)];

    theta_new=theta-alpha*(H\grad);

    %when to stop
    if norm(theta_new-theta)<=delta
        theta=theta_new;
        
        break
    end

    theta=theta_new;
end

f_hat=theta(1)
b_hat=theta(2)

    %idpoly model and validation


    A=1;
    B=[zeros(1,nk) b_hat];
    C=1;
    D=1;
    F=[1 f_hat];
    Ts=0.01;

    model_oe=idpoly(A,B,C,D,F,0,Ts);

    start_step=150+200+200+1;
    end_step=start_step+70-1;

    u_val=u(start_step:end_step);
    y_val=y_data(start_step:end_step);


    u_val=u_val/su;
    y_val=y_val/sy;

    u_val=u_val(:);
    y_val=y_val(:);
    data_val=iddata(y_val,u_val,Ts);

    figure(4)
    compare(data_val,model_oe);
