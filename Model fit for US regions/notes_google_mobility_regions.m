[a,b,c]=xlsread('retail_and_recreation_percent_change_from_baseline_ST.csv');
dates=c(1,5:end);
dates_n=datenum(dates);
states_regions=a(:,4);
states_populations=a(:,3);
states_names=c(2:end,2);
google_mobility_retail_recreation=a(:,5:end);
%%
a1=(100+google_mobility_retail_recreation)./100;  
for m=1:4;
    i1=find(states_regions==m);
    pop_tot=sum(states_populations(i1));
    aux=states_populations(i1)./pop_tot; 
    aux1=repmat(aux, 1, length(a1));
    table_regions(m,:)=sum(aux1.*a1(i1,:),1); 
end;
%%
a2=movmean(table_regions',7)'; 
figure; plot(dates_n, a2, '-');
datetick('x');
%%
% checks
% sum(isnan(google_mobility_retail_recreation(:)))
% i1=find(states_regions==4); whos i1
% states_names(i1)
%%
save retail_and_recreation_percent_change_from_baseline_ST.mat