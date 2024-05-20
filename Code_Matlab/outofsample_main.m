clear all;
clc;

addpath('functions')
addpath('data')
addpath('Historical FRED-MD Vintages Final')


%% Choices
models=[3 4 5 6 7 ];%select models 
% 3:Horseshoe 4:ridge 5:lasso 6:Gaussian Processes 7:Random forest
quantgrid= [0.05 .1  .25 .5 .75  .9 .95];% grid of quantiles

selectVariable=1;% 1:CPI 2:Indpro 3:umcs 4:employment 

fred=0;% if 1 use fred data
text=1;%0: use no text data 1:use text data with 80 topics 2:use text data with 80 topics combined with sentiment
hfore=2;% 1:nowcast 2:one step ahead forecast
ylag=12;%select number of lags for the dependend variables
xlag=1; %select the number of lags for the predictors must be >0


%%
if selectVariable==1
tcode_opt = 6; 
select_target="CPIAUCSL";%"PCEPI"
vintage = 1;
elseif selectVariable==2
select_target="INDPRO";
tcode_opt = 5; 
vintage = 1;
elseif selectVariable==3
select_target="UMCSENTx";
tcode_opt = 2;
vintage = 0;
elseif selectVariable==4
tcode_opt = 5; 
select_target="CE16OV";%"PCEPI
vintage = 1;
end

% 1.6.1980 insample starts
% 1999:08 first out of sample forecast 

startsample_month=258;
timelag=0;%do not change
%%
Input.delta=1;% 1: constant variance <1: discouting
Input.epsilon=0.0001;
Input.first=1;

%% Target var

[data1, header1]=xlsread('2022-01.csv');
[~, header3]=xlsread('1999-08.csv');
vars_lastvintage = string(header1(1,2:end));
vars_1999 = string(header3(1,2:end));
select =contains(vars_lastvintage, vars_1999);
obj_vars = vars_lastvintage(select);
targe_var_str =select_target; % Use this to select which variable or we can also specify the string by name
target_var_select = ismember(vars_lastvintage, targe_var_str);

%% datum


date_mont=char(header1(231+2-2+259:end,1)); %start in 1999:08
date_mont_char=[date_mont(:,7:end),repmat('-', size(date_mont,1),1),date_mont(:,1:2)];

%% macro data

head3=string(header3(1,2:end));
macronames=head3([1:61 64:65 88:108]);
transcode=data1(1,:);
head1=string(header1(1,2:end));
select=contains(head1,string(macronames));
macronames=head1(select);
tcode_macro=transcode(select);


%% text data

minus=0;
if text==1 %start in 06.1980 end in 12M2021 
    [data, header_text]=xlsread('ctm_80_10_topicsonly_monthly.csv');
    data_text=data(1:end-minus,:); %remove 2021
elseif text==2
     [data, header_text]=xlsread('ctm_80_10_senttopics_bignomics_monthly.csv');
    data_text=data(1:end-minus,:); %remove 2021
elseif text==3

elseif text==4

elseif text==5
     [data, header_text]=xlsread('sentiment_barbaglia.csv');
    data_text=data(1:end-minus,:); %remove 2021
 startsample_month =   startsample_month+43;
elseif text==6
        [data, ~]=xlsread('ctm_80_10_topicsonly_monthly.csv');
    data_text1=data(1:end-minus,:); %remove 2021
         [data, ~]=xlsread('ctm_80_10_senttopics_bignomics_monthly.csv');
    data_text2=data(1:end-minus,:); %remove 2021
        data_text=[ data_text1  data_text2];
        elseif text==7
        [data, ~]=xlsread('ctm_68_10_topics_monthly.csv');
    data_text1=data(1:end-minus,:); %remove 2021
         [data, ~]=xlsread('ctm_68_10_senttopics_bignomics_monthly.csv');
    data_text2=data(1:end-minus,:); %remove 2021
        data_text=[ data_text1  data_text2];
        elseif text==8
        [data, ~]=xlsread('ctm_100_10_topics_monthly.csv');
    data_text1=data(1:end-minus,:); %remove 2021
         [data, ~]=xlsread('ctm_100_10_senttopics_bignomics_monthly.csv');
    data_text2=data(1:end-minus,:); %remove 2021
        data_text=[ data_text1  data_text2];
        elseif text==9
        [data, ~]=xlsread('ctm_80_10_FullSample_topics_monthly.csv');
    data_text1=data(1:end-minus,:); %remove 2021
         [data, ~]=xlsread('ctm_80_10_FullSample_senttopics_bignomics_monthly.csv');
    data_text2=data(1:end-minus,:); %remove 2021
        data_text=[ data_text1  data_text2];
        elseif text==10
        [data, ~]=xlsread('ctm_80_10_topics_weekly.csv');
    data_text1=data(1:end-minus,:); %remove 2021
         [data, ~]=xlsread('ctm_80_10_senttopics_bignomics_weekly.csv');
    data_text2=data(1:end-minus,:); %remove 2021
        data_text=[ data_text1  data_text2];
elseif text==11
            [data, ~]=xlsread('ctm_80_10_topics_weekly.csv');
    data_text1=data(1:end-minus,4:4:80*4); %remove 2021
         [data, ~]=xlsread('ctm_80_10_senttopics_bignomics_weekly.csv');
    data_text2=data(1:end-minus,4:4:80*4); %remove 2021
        data_text=[ data_text1  data_text2];
                elseif text==12
        [data, ~]=xlsread('ctm_80_15_topics_monthly.csv');
    data_text1=data(1:end-minus,:); %remove 2021
         [data, ~]=xlsread('ctm_80_15_senttopics_bignomics_monthly.csv');
    data_text2=data(1:end-minus,:); %remove 2021
        data_text=[ data_text1  data_text2];
elseif text==13
           [data, ~]=xlsread('stm_80_10_topics_monthly.csv');
    data_text1=data(1:end-minus,:); %remove 2021
         [data, ~]=xlsread('stm_80_10_senttopics_bignomics_monthly.csv');
    data_text2=data(1:end-minus,:); %remove 2021
        data_text=[ data_text1  data_text2]; 
end






%ende M12.2020
%% fred financial data

financialnames=string(head3([62:63 66:87]));
% Exclude possible objective variable
ex_fin= ~contains(string(financialnames), targe_var_str);
select=contains(head1,financialnames(ex_fin));
financialdata=data1(startsample_month+1:end,select);%start in 07.1980 end in 12M2020
nomissings=(sum(ismissing(financialdata)))==0;
financialdata=financialdata(:,nomissings);
tcode=transcode(select);
tcode=tcode(nomissings);
financialdata_trans=transform(financialdata,tcode);
financialdata_trans=financialdata_trans(2:end,:);% %start in 08.1980 end in 12M2020 becuase of tcode=5

% dates_financialdata = original_dates(2:end,:);

%ende M12.2020 auf der vintageDatei M1.2021
%% target VAR
targetvar_levels = data1(startsample_month+1:end,target_var_select);




%% out of sample
for kk=quantgrid
    Input.quant=kk;

    for ii=models

        T=size(date_mont,1)-hfore+1;
        saveoutput=NaN(T,2);

        for tt=1:T
            disp([num2str( 100*tt/T) '% of forecasts done. Model: ' num2str(find(models == ii)) ...
                ' out of ' num2str(size(models,2)) ...
                '. Quant: ' num2str(find(quantgrid == kk)) ' out of ' num2str(size(quantgrid,2))])

            %% fred macro month
                [data2, header2] = xlsread([date_mont_char(tt,:), '.csv']);

                % Exclude the target var
                head=string(header2(1,2:end));


            if fred>0


                exclude = ~contains(string(macronames), targe_var_str);
                macronames_exclude = string(macronames);
                selectmacro=contains(head, macronames_exclude(exclude));
                if sum(selectmacro)>size(head,2)
                    error('something is wrong with selectmacro')
                end

                tcode=tcode_macro(contains(macronames_exclude,head));
                macrodata=data2(startsample_month+1:end,selectmacro);
                macrodata_trans=transform(macrodata,tcode);
                if tcode_opt==6
                macrodata_trans=macrodata_trans(3:end,:);%remove first two because of tcode=6
                else
                macrodata_trans=macrodata_trans(2:end,:);%remove first two because of tcode=6   
                end
                missings=sum(ismissing(macrodata_trans(end-4:end,:)));
                maxmissing=max(missings);
                dataNOmissings=NaN(size(macrodata_trans,1)-maxmissing,size(macrodata_trans,2));
                for iii=1:size(macrodata_trans,2)
                    dataNOmissings(:,iii)=macrodata_trans(1+maxmissing-missings(iii):end-missings(iii),iii);
                end
                dataNOmissings= dataNOmissings(:,(sum(ismissing(dataNOmissings)))==0);
                Tnomis=size(dataNOmissings,1);
                Xmacro=dataNOmissings(end-Tnomis+1:end,:);

                % combine macro and finance

                diff=(T-tt)+timelag;
                financialdata2=financialdata_trans(1:end-diff,:);
                Xfinancial= financialdata2(1:end,:);
                diff=size(Xfinancial,1)-size(Xmacro,1);
                Xfred=[Xmacro Xfinancial(1+diff:end,:)];
            else
                Xfred=[];
            end
            %% text
            if text>0
                diff=(T-tt)+timelag;
                Xtext= data_text(1:end-diff,:);
                if fred>0
                    diff=size(Xtext,1)-size(Xfred,1);
                    Xtext=Xtext(1+diff:end,:);
                end
            else
                Xtext=[];
            end

            %%


            Xmonth=[Xfred Xtext];

            if xlag>1
                Xmonthlag=mlag2(Xmonth,xlag-1);
                Xmonth=[Xmonth(xlag:end,:) Xmonthlag(xlag:end,:)];
            end

            

            %% Objective variable
            % Set the target variable
                if vintage == 1
                include = contains(head, select_target);
                y_or=data2(startsample_month+1:end,include);
                else
                y_or = targetvar_levels(1:end-T+tt-hfore);

                end



                y = yfcst2(y_or,tcode_opt,hfore,0);
                actualvalue = yfcst2(targetvar_levels(1:end),tcode_opt,hfore,0);
                actualvalue=actualvalue(end-T+tt-hfore,1);

            if ylag>0
                if tcode_opt == 5
                    yx= y_or;
                    yx = yfcst(yx,5,1); 
                    if ylag>1
                    y_lag=[yx(1:end-1) mlag2(yx(1:end-1),ylag-1)];
                    y_lag=y_lag(ylag:end,:);
                    else
                    y_lag=yx(1:end-1);
                    end
                    
                elseif tcode_opt == 6
                    yx= y_or;
                    yx = yfcst(yx,6,1);
                    if ylag>1
                    y_lag=[yx(2:end-1) mlag2(yx(2:end-1),ylag-1)]; 
                    y_lag=y_lag(ylag:end,:);
                    else
                    y_lag=yx(2:end-1);
                    end
                    
                elseif tcode_opt==2
                    yx= y_or;
                    yx = yfcst(yx,2,1);
                    if ylag>1
                    y_lag=[yx(1:end-1) mlag2(yx(1:end-1),ylag-1)];
                     y_lag=y_lag(ylag:end,:);
                     else
                    y_lag=yx(1:end-1);
                    end
                   
                elseif tcode_opt==4
                    yx= y_or;
                    yx = yfcst(yx,4,1);
                    if ylag>1
                    y_lag=[yx(1:end-1) mlag2(yx(1:end-1),ylag-1)];
                    y_lag=y_lag(ylag:end,:);
                    else
                    y_lag=yx(1:end-1);
                    end
                    
                end
                if fred>0 || text>0
                    diff=max(size(Xmonth,1)-size(y_lag,1),0);
                    Xmonth=Xmonth(1+diff:end,:);

                    diff=max(size(y_lag,1)-size(Xmonth,1),0);
                    y_lag=y_lag(1+diff:end,:);
                    X=[y_lag Xmonth];
                else
                    X=y_lag;
                end
            else
                X=Xmonth;
            end


          


            diff=max(size(y,1)-size(X,1),0);
            y=y(1+diff:end,1);
            diff=max(size(X,1)-size(y,1),0);
            X=X(1+diff:end,:);
            [zz,~,~]=transform_data(X,2);%transform to mean zero std one

            if ii==8 || ii==9
            [~,zzz,~,~,~] = factors_em(X(:,ylag:end),3,2,2);
            [zzzz,~,~]=transform_data(X(:,1:ylag),2);%transform to mean zero std one
             zz=[zzzz zzz];
            end


            if ii==6
            zz=  cholkernel(zz);
            end

            z=[ones(size(zz,1),1) zz];

            Input.zin=z(1:end-hfore,:);
            Input.zout=z(end,:);
            [Input.yt,yt_mut,yt_sdt]=transform_data(y(1:end-hfore,1),2);



            if ii==1
                [p_beta,~]=quantreg(Input.zin,Input.yt,Input.quant);
                Output.pointf=     Input.zout*p_beta;
            elseif ii==2
                Output=regress(Input);
            elseif ii==3
                Output=horseshoe(Input);
            elseif ii==4
                Output=ridge(Input);
            elseif ii==5
                Output=LASSO(Input);
            elseif ii==8
                [p_beta,~]=quantreg(Input.zin,Input.yt,Input.quant);
                Output.pointf=     Input.zout*p_beta;
            elseif ii==9
                Output=regress(Input);
            elseif ii==6
                Output=regresskernel(Input);
            elseif ii==7
                rng(1); % For reproducibility
                Mdl = TreeBagger(100,Input.zin,Input.yt,'Method','regression');
             Output.pointf = quantilePredict(Mdl,Input.zout,'Quantile',Input.quant);
            end



            
            Output.pointf=Output.pointf*yt_sdt(1)+yt_mut(1);

            saveoutput(tt,1)=actualvalue;
            saveoutput(tt,2)=(Output.pointf);

        end
if text==5
minusstart=43;
else
minusstart=0;
end
                temp=['Month_pred','Target=',targe_var_str,'quant=',num2str(Input.quant),'hfore=',num2str(hfore),...
                    'model=',num2str(ii),'fred=',...
                    num2str(fred),'text=',num2str(text), 'ylag=',num2str(ylag),...
                    'start=',num2str(startsample_month-minusstart),'tcode=',num2str(tcode_opt),'.mat'];
                temp = strjoin(temp, ' ');
                Quant=Input.quant;
                save(temp,'saveoutput','Quant');
    end
end