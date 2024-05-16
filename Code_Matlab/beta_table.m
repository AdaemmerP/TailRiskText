% Extract the betas
clear all;
clc;

addpath('functions')
addpath('data')
addpath('Historical FRED-MD Vintages Final')
addpath('Results')


%% Select models and quantile
modelnames={ 'AR(1)','AR1 Bayes','HS','RIDGE','LASSO', 'BayesianKernel', 'Tree'};
nmodels= [3 4 5 6 7 ];% [3 4];%
modelnames = modelnames(nmodels);
quantgrid= [0.1 0.5 0.90];%  0.5; %We can specify more quants
quantnames= {'10%' '50%' '90%'}; %"50%"; %

%% Basic variables
Input.epsilon=0.0001;
Input.delta=1;
Input.first=0;
addconst=1;
startsample_month=258+1-1;
start_outofsample = 490;
start_outofsample_text = 234;
hfore=2;
ylag=12;
xlag=1;
text_option = 6;


%% Select objective variables and information
select_variables =   [1 2 3 4];%5 7  % 1 ; %  We can specify one but more for the loop
tcode_options = [6 5 2 5  11 11 11];
variable_names = ["CPIAUCSL", "INDPRO",  "UMCSENTx", "CE16OV","UNRATE", "HOUST",...
     "PCEPI"];
vintage_options = [1 1 1 0 0 0 1];
tcode_options = tcode_options(select_variables);
target_var_select = variable_names(select_variables);
vintage_options = vintage_options(select_variables);
%% Load FRED data
[data, header]=xlsread('2022-01.csv');
[~, header3]=xlsread('1999-08.csv');
head3=string(header3(1,2:end));
head1=string(header(1,2:end));
transcode = data(1,:);
macronames= head3([1:61 64:65 88:108]);
financialnames=head3([62:63 66:87]);
date_mont=char(header(231+2-2+259:end,1)); %out of sample % Jan revision: -2 minus 2 out of sample
date_mont_char=[date_mont(:,7:end),repmat('-', size(date_mont,1),1),date_mont(:,1:2)];
T=size(date_mont,1)-hfore+1;
% Remove the first periods until the first out-of-sample observation
data = data(1+startsample_month:end,:);

dataselect_names = {'Fred' 'Text' 'Both'};

% Dimension ordering
% Change here for ordering: quantiles(Rows), betas, dataselect (subplot columns), models (Columns)
Beta_save=NaN(size(quantgrid, 2), 2000,  size(nmodels,2));
length_save=NaN(size(nmodels,2),1); % Models
beta_lengthsave = NaN(3,size(target_var_select,2)); % Data select, 1, variables
beta_names = strings(2000,3*size(quantgrid,2)); % Betas, Data select* quantgrid - Now one at a time


%% Names of all possible betas
% Macro names
select = contains(head1, macronames);
allmacro = head1(select);

% All financial
select = contains(head1, financialnames);
allfin = head1(select);
fin_vars = data(1:end, select);
nomissings = (sum(ismissing(fin_vars)))==0;
allfin = allfin(nomissings);
all_posvars = [allmacro  allfin];


% Check the objective variables
select = ~contains(variable_names, all_posvars);
all_posvars = [all_posvars variable_names(select)];

% Text data
minus=0;
if text_option==1 %start in 07.1980 end in 12M2021 
    [datat, header_text]=xlsread('topics_all_K80_10TsdWrds_No_compounds.csv');
    data_text=datat(1:end-minus,:); %remove 2021
elseif text_option==2
     [datat, header_text]=xlsread('topics_all_K100_10TsdWrds_No_compounds_new.csv');
    data_text=datat(1:end-minus,:); %remove 2021
elseif text_option==3
     [datat, header_text]=xlsread('topics_all_K80_15TsdWrds_No_compounds_new.csv');
    data_text=datat(1:end-minus,:); %remove 2021
elseif text_option==4
     [datat, header_text]=xlsread('topics_all_K100_15TsdWrds_No_compounds_new.csv');
    data_text=datat(1:end-minus,:); %remove 2021
elseif text_option==6
            [datat, header_text1]=xlsread('ctm_80_10_topicsonly_monthly.csv');
    data_text1=datat(1:end-minus,:); %remove 2021
         [datat, header_text2]=xlsread('ctm_80_10_senttopics_bignomics_monthly.csv');
    data_text2=datat(1:end-minus,:); %remove 2021
        data_text=[ data_text1  data_text2];
end
if text_option==6
text_names=[ string({header_text1{1,2:end}})  string({header_text2{1,2:end}})];
else
text_names = string({header_text1{1,2:end}});
error("check code")
end
% One big list of names
if xlag > 1
    all_betas=string([strcat('ylag-', string(1:ylag)) ...
        strcat(repmat(all_posvars, [1, xlag]), ...
        sort(repmat(strcat('-lag', string(1:xlag)), size(all_posvars)))) ...
        strcat(repmat(text_names, [1, xlag]), ...
        sort(repmat(strcat('-lag', string(1:xlag)), size(text_names))))]);
else
    all_betas=string([ strcat('ylag-', string(1:ylag)) ...
        strcat(all_posvars) ...
        strcat(text_names)]);
end


% Names of model combinations

comb_names = strings([(3 * size(nmodels,2)) 1]);


%%%%%%%% Loop starts for: 1. quant 2. variable 3. Type of data 4. Number of
%%%%%%%% models specified outer loops its used also for exporting the plots

progress = size(target_var_select,2) * 1 * size(nmodels,2) * size(quantgrid,2);
progress_counter = 0;

countvariable = 0;
for variable =  target_var_select
    countvariable = countvariable + 1;
    tcode_opt = tcode_options(countvariable);
    countquantile = 0;
    countbetanames = 0;
    for qq = quantgrid
        countquantile = countquantile+1;
        for dataselect = 3% 1:3 vorher
            if dataselect==1
                fred=1;
                text2=0;
            elseif dataselect==2
                fred=0;
                text2=1;
            elseif dataselect==3
                fred=1;
                text2=6;
            end
            countbetanames = countbetanames + 1;
            %% Fred data
            if fred == 1
                % Exclude objective variable
                exclude = ~contains(head1,variable);
                all_vars = head1(exclude);
                xdata = data(:,exclude);
                transcode2 = transcode(exclude);

                % Macro vars
                select = contains(all_vars, macronames);
                macrovars = xdata(:, select);
                macronames_exclude = all_vars(select);
                tcode_macro = transcode2(select);
                macrodata_trans=transform(macrovars,tcode_macro);
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
                dataNOmissings2= dataNOmissings(:,(sum(ismissing(dataNOmissings)))==0);
                macronames_exclude  = macronames_exclude(:,(sum(ismissing(dataNOmissings)))==0);
                Tnomis=size(dataNOmissings2,1);
                Xmacro=dataNOmissings2(end-Tnomis+1:end,:);

                % Financial variables
                select = contains(all_vars, financialnames);
                financial_vars = xdata(1:end, select);
                nomissings = (sum(ismissing(financial_vars)))==0;
                financial_vars = financial_vars(:,nomissings);
                head_financial = all_vars(select);
                financialnames2 = head_financial(nomissings);
                tcode_fin = transcode2(select);
                tcode_fin = tcode_fin(nomissings);
                financial_vars_trans = transform(financial_vars, tcode_fin);
                financial_vars_trans = financial_vars_trans(2:end,:); % Because of tcode = 5

                % Merge
                diff=size(financial_vars_trans,1)-size(Xmacro,1);%tcode6 but Tnomis has new changes
                % Xfred = [Xmacro financial_vars(1+diff:end,:)];
                Xfinancial = financial_vars_trans(1+diff:end,:); % new

            else
                Xmacro = [];
                Xfinancial = [];
            end

            %% Text data
            if text2 > 0
                Xtext = data_text;
                if fred>0
                    diff=size(Xtext,1)-size(Xmacro,1);
                    Xtext=Xtext(1+diff:end,:);
                end
            elseif text2 == 0
                Xtext = [];
            end
            %% Financial and text
            Xfintex = [Xfinancial Xtext];

            %% Lag X
            if xlag>1
                Xmacrolag = mlag2(Xmacro, xlag-1);
                Xmacro = [Xmacro(xlag:end,:) Xmacrolag(xlag:end,:)];
                Xfintexlag = mlag2(Xfintex, xlag-1);
                Xfintex = [Xfintex(xlag:end,:) Xfintexlag(xlag:end)];
            end

            %% Lagged objective variable
            select = contains(head1, variable);
            y_or = data(:,select);

            if ylag>0
                if tcode_opt == 5
                    yx= y_or;
                    yx = yfcst(yx,5,1);
                    y_lag=mlag2(yx(1:end),ylag); % Deleted the (1:end-1)
                    y_lag=y_lag(1+ylag:end,:);
                elseif tcode_opt == 6
                    yx= y_or;
                    yx = yfcst(yx,6,1);
                    y_lag=mlag2(yx(2:end),ylag); % Deleted the (1:end-1)
                    y_lag=y_lag(1+ylag:end,:);
                elseif tcode_opt==2
                    yx= y_or;
                    yx = yfcst(yx,2,1);
                    y_lag=mlag2(yx(1:end),ylag); % Deleted the (1:end-1)
                    y_lag=y_lag(1+ylag:end,:);
                elseif tcode_opt==4
                    yx= y_or;
                    yx = yfcst(yx,4,1);
                    y_lag=mlag2(yx(1:end),ylag); % Deleted the (1:end-1)
                    y_lag=y_lag(1+ylag:end,:);
                end

                if fred>0 || text2>0
                    diff=max(size(Xmacro,1)-size(y_lag,1),0);
                    Xmacro=Xmacro(1+diff:end,:);
                    diff=max(size(Xfintex,1)-size(y_lag,1),0);
                    Xfintex=Xfintex(1+diff:end,:);
                    diff=max(size(y_lag,1)-size(Xmacro,1),0);
                    y_lag=y_lag(1+diff:end,:);
                    diff=max(size(y_lag,1)-size(Xfintex,1),0);
                    y_lag=y_lag(1+diff:end,:);
                    X=[y_lag Xmacro];
                else
                    X=y_lag;
                end
            else
                X= Xmacro;
            end

            % Direct forecasting correction
            X = [X(1:end-hfore,:) Xfintex(2:end-hfore+1,:)]; % Correction for fin and text vars
            Input.quant=qq;

            % % Names of betas with the lags
            if dataselect ==1
                if xlag > 1
                    txt=string([ strcat('ylag-', string(1:ylag)) ...
                        strcat(repmat(macronames_exclude, [1, xlag]), ...
                        sort(repmat(strcat('-lag', string(1:xlag)), size(macronames_exclude)))) ...
                        strcat(repmat(financialnames2, [1, xlag]), ...
                        sort(repmat(strcat('-lag', string(1:xlag)), size(financialnames2))))]);
                else
                    txt=string([ strcat('ylag-', string(1:ylag)) ...
                        strcat(macronames_exclude) ...
                        strcat(financialnames2)]);
                end

            elseif dataselect==2
                if xlag>1
                    txt=string([strcat('ylag-', string(1:ylag)) ...
                        strcat(repmat(strcat({header_text{1,2:end}}), [1, xlag]), ...
                        sort(repmat(strcat('-lag', string(1:xlag)), size(strcat(text_names)))))]);
                else
                    txt=string([strcat('ylag-', string(1:ylag)) ...
                         {header_text{1,2:end}}]);
                end

            elseif dataselect==3
                if xlag > 1
                    txt=string([strcat('ylag-', string(1:ylag)) ...
                        strcat(repmat(macronames_exclude, [1, xlag]), ...
                        sort(repmat(strcat('-lag', string(1:xlag)), size(macronames_exclude)))) ...
                        strcat(repmat(financialnames2, [1, xlag]), ...
                        sort(repmat(strcat('-lag', string(1:xlag)), size(financialnames2)))) ...
                        strcat(repmat(strcat(text_names), [1, xlag]), ...
                        sort(repmat(strcat('-lag', string(1:xlag)), size(strcat(text_names)))))]);
                else
                    txt=string([ strcat('ylag-', string(1:ylag)) ...
                        strcat(macronames_exclude) ...
                        strcat(financialnames2) ...
                        text_names]);
                end
            end
            beta_names(1:size(txt,2), countbetanames) = txt';
            beta_lengthsave(dataselect,1) = size(txt,2);

            countmodel=0;
            for ii = nmodels
                countmodel=countmodel+1;
               % name=['Month_pred','Target=',variable,'quant=', ...
                %    num2str(Input.quant),'hfore=',num2str(hfore),...
                 %   'model=',num2str(ii),'fred=',...
                  %  num2str(fred),'text=',num2str(text2), 'ylag=',num2str(ylag),...
                   % 'xlag=',num2str(xlag),'tcode=',num2str(tcode_opt),'.mat'];
                name=['Month_pred','Target=',variable,'quant=',num2str(Input.quant),'hfore=',num2str(hfore),...
                    'model=',num2str(ii),'fred=',...
                    num2str(fred),'text=',num2str(text2), 'ylag=',num2str(ylag),...
                    'start=',num2str(startsample_month),'tcode=',num2str(tcode_opt),'.mat'];
                name = strjoin(name, ' ');
                load(name)

                length_save(countmodel, 1) = size(X,2); % Length number of columns
                quantforecast = saveoutput(:,2);
                [cutobs,~]=size(quantforecast);
                xx=X(end-cutobs+1:end,:);
                [z,~,~]=transform_data(xx,2);
                z=[ones(size(z,1),1) z];
                Input.zin=z;
                Input.zout=z(end,:);
                [Input.yt,yt_mut,yt_sdt]=transform_data(quantforecast,2);
                %Output=horseshoemcmc(Input);
                %bb=Output.betamean(2:end)
                % Change here for ordering: quantiles(Rows), betas, dataselect (subplot columns), models (Columns)
                rng default % For reproducibility
                [B,FitInfo] = lasso(z,Input.yt,'CV',10,'Intercept',false,'Standardize',false);
                idxLambdaMinMSE = FitInfo.IndexMinMSE;
                bb= B(2:end,idxLambdaMinMSE );
                Beta_save(countquantile, 1:size(z,2)-1, countmodel)=bb;
                progress_counter = progress_counter + 1;
                disp([num2str(progress_counter) ' out of ' num2str(progress) ...
                    ' single_var_dataselect'])

            end
        end
    end

    % Correction to get the betanames
    beta_names2 = beta_names(:,1:1);
    length_beta = NaN(size(beta_names2,2),1);
    for i = 1:size(beta_names2,2)
        length_beta(i,1) = sum(beta_names2(:,i) ~= "");
    end
    beta_names2 = beta_names2(1:max(length_beta),:);

    %%%%% New loops for the beta tables
    for qq2 = 1:size(quantgrid,2)
        big_beta_cont = NaN(size(all_betas,2),size(nmodels,2)*1);
        % Names of model combinations
        comb_names = strings([(1 * size(nmodels,2)) 1]);
        for model = 1:size(nmodels,2)
            % First column of beta names
            test = beta_names2;
            test_values = squeeze(Beta_save(qq2,1:max(length_beta),model));
            %             ind = strings([size(test,1) 3]);
            beta_ord2 = NaN([size(all_betas,2) 1]);
            for j = 1:size(test,2) % Loop columns
                for i = 1:size(test,1) % Loop over the rows
                    aux = find(matches(all_betas, test(i,j)));
                    %                     if isempty(aux)
                    %                         ind(i,j) = "";
                    %                     else
                    %                         ind(i,j) = aux;
                    %                     end
                    beta_ord2(aux,j) = test_values(i);
                end
            end
            big_beta_cont(:,(1:1)+(model-1)*1) = beta_ord2;
            comb_names((1:1)+(model-1)*1, 1) = strcat("Model=" , ...
                modelnames(model),"_",  "_", ...
                'ylag=',num2str(ylag),"_",...
                'xlag=',num2str(xlag),"_",'tcode=',num2str(tcode_opt));
        end

        bib_beta_cont = array2table(big_beta_cont', 'VariableNames', all_betas, ...
            'RowNames', comb_names);
        sheet_name = strcat(variable,"q=", quantnames(qq2), "hfore", ...
            num2str(hfore));
        writetable(bib_beta_cont,"beta_results.xlsx",'Sheet',sheet_name, ...
            'WriteRowNames',true)
    end
end


