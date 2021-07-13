[a,b,c]=xlsread('91-DIVOC-deaths-regions-normalized.csv');
deaths_regions=a;
dates_deaths=c(1,2:end);
dates_deaths_n=datenum(dates_deaths);
region_names{1}='South US'; region_names{2}='West US'; 
region_names{3}='Northeast US'; region_names{4}='Midwest US';
%%
clear fw* sw* *array*
for m=1:4
deaths_curr=deaths_regions(m,:)';
[BETA_FW,R_FW,J_FW,COVB_FW,MSE_FW] = ...
    nlinfit((1:140)',deaths_curr(1:140)/(IFR.*10^5),@(beta,X)SIR_CIRfit(beta,[0.95,1.92,230],X),[0.0001,0.95,15,1.15]);
BETA_FW_CI= nlparci(BETA_FW,R_FW,'covar',COVB_FW);
[BETA_SW,R_SW,J_SW,COVB_SW,MSE_SW] = ...
    nlinfit(([140:350])',deaths_curr([140:350])/(IFR.*10^5),@(beta,X)SIR_CIRfit(BETA_FW,beta,X),[0.95,1.92,230]);
%
BETA_SW_CI= nlparci(BETA_SW,R_SW,'covar',COVB_SW);
%
[RM_sample, inc_sample]=SIR_CIR_output(BETA_FW,BETA_SW,1:390);
%
figure;
plot(dates_deaths_n, deaths_curr, 'ko');
hold on; plot(datenum(2020,3,31)+(1:350),IFR*10^5.*inc_sample(1:350),'k-');
datetick('x');
title(region_names{m});
xlabel('Date in 2020-2021');
ylabel('Daily deaths per 100,000');
%%
RM_array(:,m)=RM_sample;
inc_array(:,m)=inc_sample;
%%
fw_params_95CI_l(1:4,m) =BETA_FW_CI(:,1);
fw_params_95CI_u(1:4,m) =BETA_FW_CI(:,2);
fw_params_bestfit(1:4,m)=BETA_FW;
%fw_params_covariance{m}=COVB_FW;
%%
sw_params_95CI_l(1:3,m) = BETA_SW_CI(:,1);
sw_params_95CI_u(1:3,m) = BETA_SW_CI(:,2);
sw_params_bestfit(1:3,m)=BETA_SW;
%sw_params_covariance{m}=COVB_SW;
end;
days_deaths_plot=datenum(2020,3,31)+(1:350);
days_RM_plot=datenum(2020,3,31)+(1:390)-25;
RM_array_1=RM_array(:,1);
RM_array_2=RM_array(:,2);
RM_array_3=RM_array(:,3);
RM_array_4=RM_array(:,4);
%
deaths_array_1=inc_array(1:350,1).*(IFR.*10^5);
deaths_array_2=inc_array(1:350,2).*(IFR.*10^5);
deaths_array_3=inc_array(1:350,3).*(IFR.*10^5);
deaths_array_4=inc_array(1:350,4).*(IFR.*10^5);