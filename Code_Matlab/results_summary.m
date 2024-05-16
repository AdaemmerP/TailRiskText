
clear;
clc;

addpath('Results')
addpath('functions')
%% choices
target=1;% 1:CPI 2:Indpro 3:umcs 4:employment 
overtime=0;% set to 1 if you want to save results over time.

nmodels=[1 4 3 5 6 7 8 2 9];%select models, note that the first model will be the benchmark model
text=[  0  2   1 6  7 8 9 11 12 13];% select text data
fred=[0 1 ];%select 1 if you want to consider fred data
Nfore=[1 2 ];%Select the forecasting horizon
lagy=[ 1  12  ];
Xlag=1 ;
Start=1;
quantgrid=[.05 .1 .25 .5 .75 .90 .95];
startsample_month=258;


savefile        = ['Results_linear_all.xlsx'];
savefile2        = ['Results_linear_Overtime.xlsx'];
%%

if target==1
 targetvar="CPIAUCSL";
tcode_opt=6;   
elseif target==2
targetvar="INDPRO";  %
tcode_opt=5;
elseif target==3
targetvar="UMCSENTx";
tcode_opt = 2;    
elseif target==4
tcode_opt = 5; 
targetvar="CE16OV";%"PCEPI"
end


comb=size(nmodels,2)*3*3;
%% works for crps only if we have all text and fred data for all models %%
for nhfore= Nfore
    for covid=1
    countquant=0;
    QL_save=NaN(size(quantgrid,2),comb,1000);
        for quant=quantgrid
            countquant=countquant+1;
           
           
           
rownames=cell(comb,1);
QLmean=zeros(comb,1);
QL_time=zeros(comb,1000);
QL_DM=zeros(comb,1);
UC_test=zeros(comb,1);
DC_test=zeros(comb,1);



count=0;
    for m=nmodels
        if m==1 ||m==2 || m==8|| m==9
            lagy=1;
        else
            lagy=12;
        end
        for Lag_y=lagy
            for Text=text
                for Fred=fred
                    for xlag=Xlag
          
      
   
             
            

count=count+1;

check=1;
try
temp=['Month_pred','Target=',targetvar,'quant=',num2str(quant),'hfore=',num2str(nhfore),...
                    'model=',num2str(m),'fred=',...
                    num2str(Fred),'text=',num2str(Text), 'ylag=',num2str(Lag_y),...
                    'start=',num2str(startsample_month),'tcode=',num2str(tcode_opt),'.mat'];
               name = strjoin(temp, ' ');
load(name)
  catch ME
    count=count-1;
    check=0;
    continue;  % Jump to next iteration of: for i
end

if check==1


if m==1
modelname=['AR1_'];
elseif m==2
modelname=['AR1 Bayes_'];
elseif m==3
modelname=['Horseshoe_'];
elseif m==4
modelname=['Ridge_'];
elseif m==5
modelname=['Lasso_'];
elseif m==6
modelname=['Gaussian Processes_'];
elseif m==7
modelname=['Random Forest_'];
elseif m==8
modelname=['factor_'];
elseif m==9
modelname=['factor Bayes_'];
else
modelname=['error']    ;
end
temp(8)=modelname;
temp(9)=[];
temp(1:7)=[];
rownames{count,1}= char(strjoin(temp, ''));

y=saveoutput(Start:end,1);
VaR=saveoutput(Start:end,2);



indalpha=y<VaR;
QL=(quant-indalpha).*(y-VaR);

if Quant~=quant
error('quant')
end
if count==1
QL_ben=QL;  
end



if covid==1
QLmean_ben=mean(QL_ben(1:end,1));

QLmean(count,1)=mean(QL(1:end,1))/QLmean_ben;

QL_time=QL_time(:,1:size(QL_ben,1));
QL_time(count,:)=cumsum(QL(1:end,1))./cumsum(QL_ben);

QL_save=QL_save(:,:,1:size(QL,1));
QL_save(countquant,count,: )=QL(1:end,1);

[~,~,~,pval_R]= test_DM_HLN(QL_ben,QL,nhfore);
QL_DM(count,1)=pval_R;

T=size(indalpha,1);
indalphasum=sum(indalpha,1);
UC=-2*log( (1-quant)^(T-indalphasum) *quant^indalphasum ) ...
    + 2*log( ( 1-indalphasum/T)^(T-indalphasum) *(indalphasum/T)^indalphasum );
UC_test(count,1)=1- chi2cdf(UC,1);

 dep=indalpha(4+1:end,1)-quant ;
 lag=[indalpha(4:end-1,1) indalpha(3:end-2,1) indalpha(2:end-3,1) indalpha(1:end-4,1)]-quant  ;
 Z= [ones(size(indalpha(1:end-4,1),1),1) lag VaR(4:end-1,1)];

 
 delta=(Z'*Z)\Z'*dep;
 DQ=delta'*(Z'*Z)*delta/(quant*(1-quant));
 DC_test(count,1)=1- chi2cdf(DQ,6);

else
    minus=24;
    QLmean_ben=mean(QL_ben(1:end-minus,1));

QLmean(count,1)=mean(QL(1:end-minus,1))/QLmean_ben;

QL_save=QL_save(:,:,1:size(QL(1:end-minus,1),1));
QL_save(countquant,count,: )=QL(1:end-minus,1);


[~,~,~,pval_R]= test_DM_HLN(QL_ben(1:end-minus,1),QL(1:end-minus,1),nhfore);
QL_DM(count,1)=pval_R;

indalpha=indalpha(1:end-minus);
T=size(indalpha,1);
indalphasum=sum(indalpha,1);
UC=-2*log( (1-quant)^(T-indalphasum) *quant^indalphasum ) ...
    + 2*log( ( 1-indalphasum/T)^(T-indalphasum) *(indalphasum/T)^indalphasum );
UC_test(count,1)=1- chi2cdf(UC,1);

 dep=indalpha(4+1:end,1)-quant ;
 lag=[indalpha(4:end-1,1) indalpha(3:end-2,1) indalpha(2:end-3,1) indalpha(1:end-4,1)]-quant  ;
 Z= [ones(size(indalpha(1:end-4,1),1),1) lag VaR(4:end-1-minus,1)];
 delta=(Z'*Z)\Z'*dep;
 DQ=delta'*(Z'*Z)*delta/(quant*(1-quant));
 DC_test(count,1)=1- chi2cdf(DQ,6);

end



end%check

                      
                  
                    end
          end
      end
   end
end

temp2=[targetvar,'q=',num2str(quant*100),'h=',num2str(nhfore)];   %,'Covid=',num2str(covid)
temp2= char((strjoin(temp2, '')));
QLmean=QLmean(1:count,:);
QL_time=QL_time(1:count,:);
QL_DM=QL_DM(1:count,:);
UC_test=UC_test(1:count,:);
DC_test=DC_test(1:count,:);


Table=table(QLmean,QL_DM,UC_test,DC_test,...
  'RowNames',rownames(1:count,:));
Table2=table(QL_time,...
  'RowNames',rownames(1:count,:));

writetable(Table,savefile,'Sheet',strcat(temp2),'WriteRowNames',true)
if overtime==1
writetable(Table2,savefile2,'Sheet',strcat(temp2),'WriteRowNames',true)
end
        end


wtails=(2* quantgrid-1).^2;
wleft=(1-quantgrid).^2;
wequal=ones(1,size(quantgrid,2));
CRPS_t=mean(2*wequal'.*QL_save(:,1:count,:),1);
CRPS_tails_t=mean(2*wtails'.*QL_save(:,1:count,:),1);
CRPS_left_t=mean(2*wleft'.*QL_save(:,1:count,:),1);

CRPS=mean(CRPS_t,3);
CRPS_tails=mean(CRPS_tails_t,3);
CRPS_left=mean(CRPS_left_t,3);
CRPS_DM=zeros(count,1);
CRPS_norm=CRPS(1,1);
CRPS_tails_norm=CRPS_tails(1,1);
CRPS_left_norm=CRPS_left(1,1);
CRPS=CRPS./CRPS_norm;
CRPS_tails=CRPS_tails./CRPS_tails_norm;
CRPS_left=CRPS_left./CRPS_left_norm;


CRPS_tails_DM=zeros(count,1);
CRPS_left_DM=zeros(count,1);

for kkk=1:count
[~,~,~,pval_R]= test_DM_HLN(squeeze(CRPS_t(1,1,:)),squeeze(CRPS_t(1,kkk,:)),nhfore);
CRPS_DM(kkk)=pval_R;
[~,~,~,pval_R]= test_DM_HLN(squeeze(CRPS_tails_t(1,1,:)),squeeze(CRPS_tails_t(1,kkk,:)),nhfore);
CRPS_tails_DM(kkk)=pval_R;
[~,~,~,pval_R]= test_DM_HLN(squeeze(CRPS_left_t(1,1,:)),squeeze(CRPS_left_t(1,kkk,:)),nhfore);
CRPS_left_DM(kkk)=pval_R;
end



CRPS=CRPS';
CRPS_tails=CRPS_tails';
CRPS_left=CRPS_left';



Table=table(CRPS,CRPS_tails,CRPS_left,CRPS_DM,CRPS_tails_DM,CRPS_left_DM,...
  'RowNames',rownames(1:count,:));
temp2=['CRPS_Series=',num2str(targetvar),'h=',num2str(nhfore),'Covid=',num2str(covid)];   
writetable(Table,savefile,'Sheet',strcat(temp2),'WriteRowNames',true)



    end
end



