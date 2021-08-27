% comparison of our mpodel with SIR 
figure; 
for m=1:4;
deaths_curr=deaths_regions(m,:)';
uc=257;
if m==4;
     uc=231;
end;
lc=140;
%%
[BETA_FW,R_FW,J_FW,COVB_FW,MSE_FW] = ...
    nlinfit((1:lc)',deaths_curr(1:lc)/(IFR.*10^5),@(beta,X)SIR_CIRfit(beta,[0.95,1.92,230],X),[0.0001,0.95,15,1.15]);
BETA_FW_CI= nlparci(BETA_FW,R_FW,'covar',COVB_FW);
%
CIR_nlm=fitnlm((lc:uc)', deaths_curr([lc:uc])/(IFR.*10^5),@(beta,X)SIR_CIRfit(BETA_FW,beta,X),[0.95,1.92,230]);
[ypred yci] = predict(CIR_nlm,(1:370)');
%%
subplot(2,2,m)
x1=[uc-1,uc,370]+datenum(2020,3,31);
y1=[0,5,5];
hold on; area(x1,y1,'FaceColor', [0.8 0.8 0.8])
%%
plot(dates_deaths_n, deaths_curr, 'kx');
plot(datenum(2020,3,31)+(1:370),IFR*10^5.*ypred,'b-');
plot(datenum(2020,3,31)+(1:370),IFR*10^5.*yci,'b--');
datetick('x');
title(region_names{m});
xlabel('Date in 2020-2021');
ylabel('Daily deaths per 100,000');
%%
%%
[BETA_FW_SIR,R_FW_SIR,J_FW_SIR,COVB_FW_SIR,MSE_FW_SIR] = ...
    nlinfit((1:lc)',deaths_curr(1:lc)/(IFR.*10^5),@(beta,X)SIR_fit(beta,[0.95,1.92,230],X),[0.0001,0.95,15,1.15]);
BETA_FW_SIR_CI= nlparci(BETA_FW_SIR,R_FW_SIR,'covar',COVB_FW_SIR);
%
SIR_nlm=fitnlm((lc:uc)', deaths_curr([lc:uc])/(IFR.*10^5),@(beta,X)SIR_fit(BETA_FW_SIR,beta,X),[0.95,1.92,230]);
[ypred_SIR yci_SIR] = predict(SIR_nlm,(1:370)');
%
plot(datenum(2020,3,31)+(1:370),IFR*10^5.*ypred_SIR,'r-');
plot(datenum(2020,3,31)+(1:370),IFR*10^5.*yci_SIR,'r--');
legend('','data','Stochastic Social Activity model','','', 'SIR model','','');
end;