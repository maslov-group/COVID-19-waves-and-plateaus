function [RM_sample, inc_sample]=SIR_output(fw_params,sw_params,time)
% The function returns two vectors: R_0*M(t) and incidence rate per capita.
% The first 4 parameters (j,R1,t1,R2) are used to fit the first wave of the 
% epidemic, prior to July 15, 2020. 
% The purpose of that procedure to it ensure that the initial
% conditions, (J,S,h) on June 1 2020 are consistent with the prior epidemic dynamics. 
% Only the last 3 parameters (R2p5,R3,t_fall,tau_s) are used for the fit 
% of the epidemic dynamics within the time-frame of interest,
% July 16, 2020-Feb 13, 2021.
% To plot Death vs. time:
%IFR=0.005;
%plot(datenum(2020,3,31)+t,IFR*10^5.*I)
%%
j=fw_params(1); 
R1=fw_params(2);
t1=fw_params(3);
R2=fw_params(4); 
%R2p5=fw_params(5); 
%
R2p5=sw_params(1); 
R3=sw_params(2);
t_fall=sw_params(3);
% R3=sw_params(1);
% t_fall=sw_params(2);
%tau_s=sw_params(4);
%
D_fall=30;
R0=2.5;
t_start=0;
time=round(time);
k=200;
dt=1./k;
Tmax=max(time)+t_start;
tau_g=5;
tau_b=365*1;
%tau_b=365*5;
%
% varied in sensitivity analysis from: 
% 18.5651 (the best fit for Northeast) 
% to 54.6216 (the best fit for West)
tau_s=30;
%tau_s=20;
%tau_s=55;
%
rw_rate=1/tau_s;
k0=0.4;
%
mu=1/(1+tau_g/tau_s);
kappa=2;
lambda=1;
%
D1=30;
t2=90; 
t2p5=120;
%t2=90+t_start;
D2=10;
D2p5=10;
t3=t_fall;
D3=D_fall;
%
h=0;
%J=k0*rw_rate.*0.001;
S=1;
att=0;
%j=0.0001;
for i=1:Tmax/abs(dt)
    M_t=1+(R1-R0)/R0.*(0.5+0.5*tanh((i*dt-t1)./D1)).*(0.5+0.5*tanh((t2-i*dt)./D2))+(R2-R0)/R0.*(0.5+0.5*tanh((i*dt-t2)./D2)).*(0.5+0.5*tanh((t2p5-i*dt)./D2p5))+(R2p5-R0)/R0.*(0.5+0.5*tanh((i*dt-t2p5)./D2p5)).*(0.5+0.5*tanh((t3-i*dt)./D3))+(R3-R0)/R0.*(0.5+0.5*tanh((i*dt-t3)./D3));
    R=M_t.*R0.*S^lambda;
    inc=M_t*S^(1+1/kappa).*j./(1+h);
dh=(M_t*j./k0-rw_rate.*h.*(1+h)).*dt;
dS=(-inc+(1-S)./tau_b)*dt;
dj=1/tau_g*((R./(1+h).^(2))-1).*j.*dt;
%h=h+dh;
h=0;
S=S+dS;
j=j+dj;
j=j*(j>0.0000001);
att=att+inc*dt;
inc=S^(1+1/kappa).*j/(1+h);
t(i)=i*dt;
inc_t(i)=inc;
RM_t(i)=M_t*R0;
end
inc_sample=inc_t(k.*(time+t_start));
inc_sample=inc_sample';
RM_sample=RM_t(k.*(time+t_start));
RM_sample=RM_sample';
end
