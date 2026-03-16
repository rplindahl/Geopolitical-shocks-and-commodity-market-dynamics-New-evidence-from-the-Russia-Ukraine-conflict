%% The Multifaceted Impact of US Trade Policy on Financial Markets

% Main Code

% DIW Berlin 21.03.2021 Lukas Boer

% Bloomberg data needs to be downloaded individually and read in
% accordingly

% For the main model the  following series are needed:
% 12- Month Treasury Yield: GB12 Govt
% 10- Year Treasury Yield: GT10 Govt
% Bloomberg Dollar Spot Index: BBDXY Index
% Russel 2000: RTY Index

% The readme file lists the tickers for the disaggregated analysis part of 
% the paper


% The main model uses "ann_base" as the event-date input "ann_chosen".
% Figure 7 can be plotted when using "ann_china" as the event-date input "ann_chosen".
% "nrep" can be lowered to use fewer draws for the bootstrap


% *************************************************************************

clear;
close all

grey = [0.5,0.5,0.5];
font_size = 8;

addpath('Geopolitical_risk_agriculture_2023/functions')

% ----------------------------------------------------------------------%
%% Read in data  
% B3588:EA3718 -> For October 2023 - March 2024

data         = xlsread('Data_food_rev.csv', 'Data_food_rev', 'B3134:EA3718');    % 04.01.2022- 01.04.2024 B3134:EA3718
data         = fillmissing(data,'previous');   % replace missing values with last observation

%%
ger10ydiff = [0; diff(data(:,14))];
%%

day          = data(:,1);   %% 2:end if variables are diffed
month        = data(:,2);
year         = data(:,3);
day_of_week  = data(:,4);               
date         = datetime(year,month,day);

cornf        = log(data(:,5));
wheatf       = log(data(:,6));
sunoil_f     = log(data(:,27));
wheats        = log(data(:,28));
corns        = log(data(:,29));
soybean_f   = log(data(:, 86));
soybean_s   = log(data(:, 87));
rice_f   = log(data(:, 88));
rapeseed_f   = log(data(:, 89));
barley_s   = log(data(:, 90));
cocoa_f   = log(data(:, 91));
coffee_f   = log(data(:, 92));
sugar_f   = log(data(:, 93));
oil_f   = log(data(:, 94));
natural_gas_f = log(data(:, 95));

natgas_ttf_m = log(data(:, 96));
%natgas_ttf_y = log(data(:, 97));
%natgas_nbp_m = log(data(:, 98));
%natgas_nbp_y = log(data(:, 99));
%natgas_psv_m = log(data(:, 100));
%natgas_psv_y = log(data(:, 101));
natgas_india = log(data(:, 102));
natgas_japan = log(data(:, 103));

%lnwheatf    = data(:,46);
%lnwheats    = data(:,47);
%lncornf    = data(:,48);
%lncorns    = data(:,49);
%lnsunf     = data(:, 50);
%lnindexown    = data(:,51);
%lneuronext    = data(:,52);
%lnsp500    = data(:,53);


%ch_wheatf   = log(data(:,25));
%ch_cornf    = log(data(:,26));
index_own    = log(data(:,24));

% Yields
euronext     = log(data(:,8));
sp500        = log(data(:,9));
ger1m        = data(:,10);
ger3m        = data(:,11);
ger1y        = data(:,12);
ger2y        = data(:,13);
ger10y        = data(:,14);
%ger30y        = data(:,15);
%fra3m        = data(:,32);
%fra1y        = data(:,34);
%fra10y        = data(:,37);
%uk3m        = data(:,40);
%uk1y        = data(:,41);
%uk10y        = data(:,44);

% Currency
useuro       = log(data(:,16));
%usdeurodif   = data(:, 64);
%ukeuro       = log(data(:,17));
%useurotest = data(:,16);

% Index
spagri       = log(data(:,21));
spenergy     = log(data(:,22));
msciworldenergy = log(data(:,23));
vstoxx       = log(data(:,30));

% Diffed
wheatfdif = diff(wheatf);
cornfdif = diff(cornf);
sunoil_fdif = diff(sunoil_f);
soybean_fdif = diff(soybean_f);
rice_fdif = diff(rice_f);
rapeseed_fdif = diff(rapeseed_f);
cocoa_fdif = diff(cocoa_f);
coffee_fdif = diff(coffee_f);
sugar_fdif = diff(sugar_f);
index_owndif = diff(index_own);
spagridif = diff(spagri);
oil_fdif = diff(oil_f);
natural_gas_fdif = diff(natural_gas_f);
natgas_ttf_mdif = diff(natgas_ttf_m);
vstoxxdif = diff(vstoxx);
euronextdif = diff(euronext);
useurodif = diff(useuro);
ger10ydif = diff(ger10y);
natgas_indiadif = diff(natgas_india);
natgas_japandif = diff(natgas_japan);

% Event dates
ann_base         = data(:,7);        % baseline events --USE--
ann_less         = data(:,18);       % Färre, största events -USE--
% ann_trade        = data(:,19);       % Endast trade relaterade
ann_attack       = data(:,20);       % Endast attacker --USE--
% ann_big          = data(:, 85);      % BIG -
ev_gpr_index    = data(:, 104);       % Test GPR index


%% Chose specification
ann_chosen     = ann_base;   % Chose event selection. _base, _less (BIG), _attack
nrep           = 50;         % bootstrap draws  
RWBS           = 0;          % 1 for Residual-wild
MBBS           = 1;          % 1 for Moving Block Bootstrap
bias_adjust_bs = 0;          % bias-adjusted via bootstrap
exog           = 0;          % add day-of-week dummies (as exog. regressors) for the reduced-form VAR estimation

c1             = 95;         % confidence bands for IRFs, 95 for 90%, 97.5 for 95%
c2             = 5;          % 5 for 90 %, 2.5 for 95%    


options        = optimoptions('fminunc','Display','none', 'TolFun',1e-30, 'TolX', 1e-30, 'MaxFunEvals', 10000, 'MaxIter', 10000);
options_fms    = optimset('Display','none', 'TolFun',1e-30, 'TolX', 1e-30, 'MaxFunEvals', 10000, 'MaxIter', 10000);
options_con    = optimoptions('fmincon','Display','none', 'TolFun',1e-30, 'TolX', 1e-30, 'MaxFunEvals', 10000, 'MaxIter', 10000);


%% Estimate 6 variable model partially identified for b1 
% as in Wright (2012)

y = [wheatf,vstoxx, euronext, ger10ydiff, useuro];  
invert_variable = 2;

[T,K]      = size(y);
p          = 5;            % Välj antal laggar

length_bl  = round(5*T^(0.25));     % block length for Bootstrap

ann_use  = ann_chosen(p+1:end,:)'; 


% estimate reduced-form VAR
[A,Sigma_hat,Uhat,V,~] = olsvarc(y,p); 
Sigma_hat              = Sigma_hat(1:K,1:K);  % Kovariansmatrisen för residualerna
Uhat                   = Uhat(1:K,:);         % 

%%
% estimate VAR w/ day-of-week dummies as exog. regressors
% 4 dummies for tuesday - friday, monday via intercept
if exog == 1
    tue = day_of_week==2; wed = day_of_week==3; thu = day_of_week==4; fri = day_of_week == 5;
    dummies = double([tue wed thu fri]);   

    Mdl            = varm(K,p);
    [VAR,~,~,Uhat] = estimate(Mdl,y, 'X', dummies);
    Uhat           = Uhat';
    results        = summarize(VAR); results.AIC
    A              = [cat(2,VAR.AR{:}); eye(K*(p-1)), zeros(K*(p-1),K)];
    V              = [VAR.Constant; zeros(K*(p-1),1)];
    Sigma_hat      = VAR.Covariance; 
    VAR.Beta;                                % dummy coefficients
end

if bias_adjust_bs == 1    % bootstrap bias-adjustment as in Wright (2012)
    rng('default');
    w = y(p+1:T,:);
    q = ones(T-p,1);     
    for j = 1:p
        q = [q y(p+1-j:T-j,:)]; 
    end
    [ahat,aadj,bias] = varest_k(y,p);  
    A_ba = aadj(2:end,:)';
    V_ba = aadj(1,:)';
    Uhat = w-(q*aadj);
    Uhat = Uhat';
    A    = [A_ba; eye(K*(p-1)), zeros(K*(p-1),K)];
    V    = [V_ba; zeros(K*(p-1),1)];
end


% Create the two distinct var-cov matrices
Uhat_t     = Uhat';          % Uhat = residuals. 
Uhat1_t    = Uhat_t(ann_use==1,:);  % Plocka ut res på eventdag
Uhat0_t    = Uhat_t(ann_use==0,:);  % Icke eventdagar
Uhat0      = Uhat0_t';        % Transpose icke event
Uhat1      = Uhat1_t';        % transpose event

T0         = length(Uhat0);      % size of event and non event              
T1         = length(Uhat1); 




Sigma1_hat = (Uhat1*Uhat1')/(T1);    % Cov för event
Sigma0_hat = (Uhat0*Uhat0')/(T0);    % Cov för icke eventsize

eigen_base = eig(Sigma1_hat/(Sigma0_hat));
eigen_base = sort(eigen_base,'descend');


[B_compl,Lambda] = SimDiag(Sigma0_hat,Sigma1_hat, 1); % calculate B and Lambdas


%Lambdas_1 = diag(Lambda);
%if B_compl(2,1) < 0
%    B_compl(:,1) = B_compl(:,1)*-1; 
%end
%if B_compl(4,2) > 0
%B_compl(:,2) = B_compl(:,2)*-1;  % adjust second shock such that negative effect on S&P47 China trade index
%end

%%

% GMM estimationx0_2
S          = vech(Sigma1_hat - Sigma0_hat);
[V0, V1,~] = cov_sigmaU(Uhat,Uhat0, Uhat1, K);    % Weighting Matrix elements, cov(Sigma_i_hat)

fun         = @(x)(S - vech(x*x'))'/(V0+V1)*(S - vech(x*x'));
x0          = ones(1,K)'/(K*10);
x0_2        = [0.02 0.01 -0.02 0 -0.02]';  % Den vi har kört


% unconstrained minimization
[x, fval]   = fminsearch(fun,x0_2,options_fms);
[x2, fval2] = fminunc(fun,x0_2, options);
[x3, fval3] = fminsearch(fun,x0,options_fms);
[x4, fval4] = fminunc(fun,x0, options);

% via constrained minimization
A_con       = [];
b_con       = [];
Aeq         = [];
beq         = [];
nonlcon     = [];       
lb          = [-Inf; -Inf; -Inf; -Inf; -Inf;; 1];
ub          = [ inf; Inf; Inf; Inf; Inf; 1];
init_b12    = x0_2;
init_c_w    = 1; 
z0          = [init_b12; init_c_w];   
     
obj         = @(x) het_obj(x,S,(V0+V1),K);
[x_c, fval_c, ~, ~] = fmincon(obj,z0, A_con,b_con,Aeq,beq,lb,ub,nonlcon,options_con);
    
x_all        = [x x2 x3 x4 x_c(1:K)];
fvals        = [fval fval2 fval3 fval4 fval_c];
[~,min_fval] = min(fvals);
fval         = min(fvals);
x            = x_all(:,min_fval);

%%
if x(invert_variable) < 0
    x = x*-1;
end

b1          = x;
normalizer  = b1(2);
b1_norm     = x/x(2)*normalizer;   
scaling_factor = 1/x(2)*normalizer;
    
B = [b1_norm(1:K), zeros(size(b1_norm(1:K),1),size(b1_norm(1:K),1)-1)];

% IRFs
h       = 100;            
horizon = 0:h;

IRF_w   = irfvar(A,B,p,h); 
IRFvec  = vec(IRF_w(1:K,:))';

IRF_compl_w   = irfvar(A,B_compl,p,h);     % using complete identification
IRFvec_compl  = vec(IRF_compl_w(1:end,:))';

% Wald Specification Tests
ann_use_m = repmat(ann_use,K,1);
Wald1 = test_sigma(Uhat,ann_use_m, T1, K);   % S'/V_w*S where V_w = (V/T1 + V/T0), assume same cov-matrix for both regimes
Wald2 = fval;                                % (S-vech(b1*b1'))'/(V0+V1)*(S-vech(b1*b1')) but with unnormed b1

%% Figure 1 Shock Series
% structural shock series
What = (b1'/Sigma_hat*Uhat)/(b1'/Sigma_hat*b1);       % Stock and Watson (2018, EJ)

days = datenum(date(1+p:end,:));
idx_ann = find(ann_chosen==1);
idx_ann_date = date(ann_chosen==1);


[tsm,tm] = downsample_ts(What,days,'monthly','sum');  
tsmn     = normalize(tsm);

%%
fig1 = figure(1);
set(fig1,'Position',[150 200 900 450]); % right pos on screen, down pos on screen, width, height
plot(tm, tsmn, '-k');
set(gca, 'YGrid', 'off', 'XGrid', 'off',  'box', 'off','FontSize',12);
datetick('x','mmmyy');
xticks([tm]);
xticklabels({'Jan17', '', '', '', '', '', 'Jul17','', '', '', '', '', 'Jan18','', '', '', '', '', 'Jul18','', '', '', '', '', 'Jan19','', '', '', '', '', 'Jul19','', '', '', '', '', 'Jan20',});
yticks([-3 -2 -1 0 1 2 3 4]);
ylim([-3 3]);
xlim([datenum(tm(1)) tm(end)]);
yl = yline(0,'-','LineWidth',0.3, 'Color', [0.6 0.6 0.6]);


% Figure Appendix A1
% daily trade policy shock series
days2 = datenum(date(1:end,:));
fig22 = figure(22);
set(fig22,'Position',[150 200 900 450]);
plot(days2(6:end), What, '-k');
set(gca, 'YGrid', 'off', 'XGrid', 'off',  'box', 'off','FontSize',9);
datetick('x','ddmmmyy');
xlim([days2(6) days2(end)]);

%% CIs for 6-var Wright Model
% MB Bootstrap for the structural impulse responses
                         
if size(y,1) == T || size(y,1) == T-1
        y_t = y';
end                        
Y = y_t(:,p:end);	
for i=1:p-1
 	Y=[Y; y_t(:,p-i:T-i)];		
end


% create Matrix of residual and proxy blocks
n_s     = ceil((T0+T1)/length_bl);       % numbers of blocks needed
U_inst  = [Uhat; ann_chosen(1+p:end)'];        % stack residuals and instrument


mbb_mat = [];                            % stack blocks
for elem = 1:(T0+T1-length_bl+1)
    step_a = U_inst(:,elem:elem+length_bl-1);
    mbb_mat = cat(3,mbb_mat,step_a);
end


IRFmat    = zeros(nrep,K*(h+1));   
IRFmat2   = zeros(nrep,K*(h+1));
IRFmat_compl = zeros(nrep,K^2*(h+1));
Wald1_r   = zeros(nrep,1);
Wald2_r   = zeros(nrep,1);
T1_num    = zeros(nrep,1);         
eig_mat   = zeros(K,nrep);         
B_compl_mat = zeros(K,K,nrep);
B_compl_r = zeros(K,K);
Lambda_mat = zeros(K,nrep);
count = 0;

 
% loop to create bootstrap samples
rng('default');
for s = 1:nrep
    s
    
    Yr       = zeros(K*p,T-p+1);
      
    % ------------------------------------------------------------------ %
    % Residual-Wild Bootstrap
    
    if RWBS == 1
    U_rc      = zeros(K,T0+T1+1);
    eta = randn(1,size(Uhat,2));  eta=repmat(eta,K,1);
	U_rc(1:K,2:T-p+1) = Uhat(1:K,:).*eta;                    % start at element 2 and use U_rc(:,2) for Y(2), Y(1) created without Residual
    T1_n     = T1;
    T0_n     = T0 + T1 - T1_n;
    
    ann_use_r = ann_use;
    
    elseif MBBS == 1
        
    % ------------------------------------------------------------------ %
    % Moving Block Bootstrap
    
    U_sample = datasample(mbb_mat,n_s,3);                  % randomly choose n_s blocks
    U_sample = reshape(U_sample,K+1,size(U_sample,2)*n_s); % reshape into one matrix
    U_sample = U_sample(:,1:T0+T1);                        % keep first T observations
    
    U_r = U_sample(1:K,:);
    U_rc = zeros(K,T0+T1);
    
    % Recenter residuals u_jl+1 = u_jl* - 1/(T-l+1)* Sum_r=0^(T-l) u_i+r
    for j = 0:n_s-1
        for i = 1:length_bl
            try
        U_rc(:,j*length_bl+i) = U_r(:,j*length_bl+i) - 1/(T0+T1-length_bl+1)*sum(Uhat(1:K,i:i+T0+T1-length_bl),2);    % exp term also via mean(mbb_mat,3)
            catch
                % continues when error: last loops cannot be run for j*length +i > T0+T1
            end
        end
    end

    T1_n     = nnz(U_sample(K+1,:));         % count non-zero elements in proxy
    
    
    % check if T1_n > 0, if == 0, draw again
    
    while T1_n == 0
        
    U_sample = datasample(mbb_mat,n_s,3);                  % randomly choose n_s blocks
    U_sample = reshape(U_sample,K+1,size(U_sample,2)*n_s); % reshape into one matrix
    U_sample = U_sample(:,1:T0+T1);                        % keep first T observations
    
    U_r = U_sample(1:K,:);
    U_rc = zeros(K,T0+T1);
    
    % Recenter residuals u_jl+1 = u_jl* - 1/(T-l+1)* Sum_r=0^(T-l) u_i+r
    for j = 0:n_s-1
        for i = 1:length_bl
            try
        U_rc(:,j*length_bl+i) = U_r(:,j*length_bl+i) - 1/(T0+T1-length_bl+1)*sum(Uhat(1:K,i:i+T0+T1-length_bl),2);    % exp term also via mean(mbb_mat,3)
            catch
                % continues when error: last loops cannot be run for j*length +i > T0+T1
            end
        end
    end

    T1_n     = nnz(U_sample(K+1,:));
    
    end
    
    U_rc = [zeros(K,1) U_rc];   % U_rc(:,1) not used as Y=A*Y(-1)+U_rc starts at Y(2)
    
    T1_num(s,1) = T1_n;
    
    T0_n     = T0+T1-T1_n;
    
    
    ann_use_r  = zeros(T1_n+T0_n,1)';        % create matrix with novel event date indices      
    ann_use_r(U_sample(K+1,:)~=0) = 1;
                             
    end                         % end of bootstrap choice
    % ------------------------------------------------------------------ %
       
    pos     = fix(rand(1,1)*(T-p+1))+1;        
	Yr(:,1) = Y(:,pos);                       % initial value for Yr chosen
    U_rc = [U_rc; zeros(K*(p-1),T-p+1)];      % companion form
       
    for i = 2:T-p+1		
        if exog == 1
            dummies_r = [dummies(pos:end,:); dummies(1:pos-p,:)];
            Yr(:,i)   = V + A*Yr(:,i-1) + U_rc(:,i) + [VAR.Beta*dummies_r(i,:)'; zeros(K*(p-1),1)]; 
        else
            Yr(:,i) = V + A*Yr(:,i-1) + U_rc(:,i);  
        end
    end
     
    yr = Yr(1:K,:);
    for i = 2:p
		yr = [Yr((i-1)*K+1:i*K,1) yr];  
    end  
    yr = yr';
        
    [Ar,Sigma_hat_r,Uhatr] = olsvarc(yr,p);        
    Uhatr       = Uhatr(1:K,:);

    if exog == 1
       Mdl             = varm(K,p);
       [VAR,~,~,Uhatr] = estimate(Mdl,yr, 'X', dummies);
       Uhatr           = Uhatr';
       Ar              = [cat(2,VAR.AR{:}); eye(K*(p-1)), zeros(K*(p-1),K)];
       Sigma_hat_r     = VAR.Covariance;  
    end    

    if bias_adjust_bs == 1    
        yboot = yr(p+1:T,:);
        qboot = ones(T-p,1);     
        for j = 1:p
            qboot =[qboot yr(p+1-j:T-j,:)]; 
        end        
        aboot = adjfunc(inv(qboot'*qboot)*qboot'*yboot,bias,K,p);      
        Uhatr = yboot-(qboot*aboot);
        Uhatr = Uhatr';
        Ar_ba = aboot(2:end,:)';
        Ar    = [Ar_ba; eye(K*(p-1)) zeros(K*(p-1),K)];
    end
        
    % ---------------------------------------- %
    % Complete identification       

    Uhatr_t  = Uhatr';
    Uhatr1_t = Uhatr_t(ann_use_r'==1,:);  
    Uhatr0_t = Uhatr_t(ann_use_r'==0,:);
    Uhatr0   = Uhatr0_t';
    Uhatr1   = Uhatr1_t';

    T0r = size(Uhatr0,2);                    
    T1r = size(Uhatr1,2);

    Sigmar1_hat = (Uhatr1*Uhatr1')/(T1r);
    Sigmar0_hat = (Uhatr0*Uhatr0')/(T0r);
    
    [B_compl_r,Lambda_bs] = SimDiag(Sigmar0_hat,Sigmar1_hat, 1);
    if B_compl_r(2,1) < 0   % For stocks
           B_compl_r(:,1) = B_compl_r(:,1)*-1;
    end
    
    Lambda_mat(:,s) = diag(Lambda_bs);

    %if B_compl_r(4,2) > 0   
    %       B_compl_r(:,2) = B_compl_r(:,2)*-1;
    %end
    B_compl_mat(:,:,s) = B_compl_r;
    % -------------------- %

    % minimization
    [z,~,z_unscaled,fval,V_inv,S_r,V0_r,V1_r] = ident_het_fminunc(ann_use_r',T0_n+T1_n,K,Uhatr,b1,options,normalizer,2,1);
    [zb,~,z_unscaledb,fvalb,~,~,~,~]          = ident_het_fmins(ann_use_r',T0_n+T1_n,K,Uhatr,b1,options_fms,normalizer,2,1);
    [z2,~,z_unscaled2,fval2,~,~,~,~]          = ident_het_fmins(ann_use_r',T0_n+T1_n,K,Uhatr,x0_2,options_fms,normalizer,2,1);
    [z2b,~,z_unscaled2b,fval2b,~,~,~,~]       = ident_het_fminunc(ann_use_r',T0_n+T1_n,K,Uhatr,x0_2,options,normalizer,2,1);
    [z3,~,z_unscaled3,fval3,~,~,~,~]          = ident_het_fmins(ann_use_r',T0_n+T1_n,K,Uhatr,x0,options_fms,normalizer,2,1);
    [z3b,~,z_unscaled3b,fval3b,~,~,~,~]       = ident_het_fminunc(ann_use_r',T0_n+T1_n,K,Uhatr,x0,options,normalizer,2,1);
    [z_unscaled4,~,fval4,~,~,~,~, eig_values] = ident_het(ann_use_r',T0_n+T1_n,K,Uhatr,z0,A_con,b_con,Aeq,beq,lb,ub,nonlcon,options_con,1);
    norm = z_unscaled4(2)/normalizer;          
    z4 = z_unscaled4/norm;   
        
    z_all_unsc = [z_unscaled z_unscaledb z_unscaled2 z_unscaled2b z_unscaled3 z_unscaled3b z_unscaled4(1:K)];
    z_all      = [z zb z2 z2b z3 z3b z4(1:K)];
    fvals = [fval fvalb fval2 fval2b fval3 fval3b fval4];
    
    count = count + (round(sum(fvals/fvals(1))/7,3)~=1.000);  
    
    [~,min_fval] = min(fvals);
    fval = min(fvals);
    z_unscaled = z_all_unsc(:,min_fval);
    z          = z_all(:,min_fval);

    if z_unscaled(invert_variable) < 0   
        z_unscaled = z_unscaled*-1;
    end
    
    b1_r    = z_unscaled;
    b1_r2   = z;

    Br  = [b1_r, zeros(size(b1_r,1),size(b1_r,1)-1)];
    Br2 = [b1_r2, zeros(size(b1_r2,1),size(b1_r2,1)-1)];
    
    % Wald Specification Tests
    ann_alt    = repmat(ann_use_r(randperm(length(ann_use_r))),K,1);   
    Wald1_r(s) = test_sigma(Uhatr,ann_alt, T1r, K);  
    Wald2_r(s) = fval;                              
    
    % eigen values from Sigma1/Sigma0
    eig_values   = sort(eig_values,'descend');
    eig_mat(:,s) = eig_values;
    
    % Compute IRFs
    IRFr           = irfvar(Ar,Br,p,h);      
    IRFr_compl     = irfvar(Ar,B_compl_r,p,h);        
    IRFmat(s,:)    = vec(IRFr(1:K,:))'-IRFvec;      
    IRFmat_compl(s,:) = vec(IRFr_compl(1:end,:))'-IRFvec_compl;   
    
    % Compute IRFs for normalized CIs
    IRFr2          = irfvar(Ar,Br2,p,h);      
    IRFmat2(s,:)   = vec(IRFr2(1:K,:))'-IRFvec;  
    
end

%%
%% -------------------------------------------------------------------- 
% Figures

% Confidence Intervals
CILO = IRFvec'-prctile(IRFmat,[c1])'; 
CIUP = IRFvec'-prctile(IRFmat,[c2])';
CILO = reshape(CILO,K,h+1); 
CIUP = reshape(CIUP,K,h+1);
CILO_compl = IRFvec_compl'-prctile(IRFmat_compl,[c1])'; 
CIUP_compl = IRFvec_compl'-prctile(IRFmat_compl,[c2])';
CILO_compl = reshape(CILO_compl,K*K,h+1); 
CIUP_compl = reshape(CIUP_compl,K*K,h+1);
CILO_2(1:K,1:h+1) = CILO(1:K,1:h+1)*100;   
CIUP_2(1:K,1:h+1) = CIUP(1:K,1:h+1)*100;
CILO_compl_2(1:K*K,1:h+1) = CILO_compl(1:end,1:h+1)*100;   
CIUP_compl_2(1:K*K,1:h+1) = CIUP_compl(1:end,1:h+1)*100;

IRF_w_2       = IRF_w(1:K,1:h+1)*100;           
IRF_compl_w_2 = IRF_compl_w(1:end,1:h+1)*100;  % complete identification, first two shocks


% Figure 2 in paper / Figure 7 if "ann_chosen = ann_china" 

name = {'Wheat future', 'VSTOXX', 'Euronext', 'Germany 10-year', 'Euro/USD Spot'}';   
labelIRF = {'%','%', '%' , ...
     '%','%', '%'}';

plot_horizons = 40;

figure(2)
for i = 1:K
    subplot(3,2,i)
    
    % Check if current subplot is for Germany 10-year
    if i == 4
        % Cumulative sums for the IRF and confidence intervals
        cumulative_IRF = cumsum(IRF_w_2(i, 1:plot_horizons));
        cumulative_CILO = cumsum(CILO_2(i, 1:plot_horizons));
        cumulative_CIUP = cumsum(CIUP_2(i, 1:plot_horizons));
        
        fill([horizon(1:plot_horizons), fliplr(horizon(1:plot_horizons))], ...
             [cumulative_CIUP, fliplr(cumulative_CILO)], ...
             'b', 'EdgeColor', 'k', 'LineStyle', '--', 'FaceAlpha', 0.1), hold on
        
        a = plot(horizon(1:plot_horizons), cumulative_IRF, 'r-', 'Linewidth', 2); hold on
    else
        % Non-cumulative plots for other subplots
        fill([horizon(1:plot_horizons), fliplr(horizon(1:plot_horizons))], ...
             [CIUP_2(i, 1:plot_horizons), fliplr(CILO_2(i, 1:plot_horizons))], ...
             'b', 'EdgeColor', 'k', 'LineStyle', '--', 'FaceAlpha', 0.1), hold on
        
        a = plot(horizon(1:plot_horizons), IRF_w_2(i, 1:plot_horizons), 'r-', 'Linewidth', 2); hold on
    end
    plot(horizon(1:plot_horizons), zeros(plot_horizons, 1), 'k', 'LineWidth', 2)
    
    title_i = strcat(name(i));
    title(title_i, 'FontSize', font_size)
    ylabel(labelIRF(i), 'FontSize', font_size)
    axis tight
    
    % Set X-axis ticks and grid
    xt = 0:5:max(horizon(1:plot_horizons));
    set(gca, 'XTick', xt);
    grid on; % This will enable both x and y grid lines
    
    % Change grid color for older MATLAB versions
    ax = gca; % get current axis handle
    set(ax, 'GridColor', [0 0 0]); % Set to a shade of gray (you can change this to your preferred color)
    
    % Set Y-axis ticks if needed (every point)
    yt = floor(min(CILO_2(i, :))):1:ceil(max(CIUP_2(i, :)));
    set(gca, 'YTick', yt);

    % Set box around the subplot
    set(ax, 'box', 'on');
    
    ax = gca; % get current axis handle
    ax.YAxis.TickLabelGapMultiplier = 0.3; % adjust the gap between y-axis and its labels
    set(gca, 'FontSize', font_size)

    title_i = strcat(name(i));
    t = title(title_i, 'FontSize', 22); % Store title handle in variable 't'
    
    ax = gca; % get current axis handle
end


% -----------------
% Figure 3 in paper
s%hift1 = sqrt(Lambda(1,1));shift2 = sqrt(Lambda(2,2));

%figure(3)
%for i = 1:K
  %  subplot(3,2,i)
   
    % Wright Impulse Responses
   % fill([horizon,fliplr(horizon)],[CIUP_2(i,1:h+1), fliplr(CILO_2(i,1:h+1))], ...
    %                    grey,'EdgeColor','none', 'FaceAlpha', 0.5), hold on
 %   a = plot(horizon, IRF_w_2(i,1:h+1),'k-', 'Linewidth', 2);  hold on
 %   b = plot(horizon, IRF_compl_w_2(i,1:h+1)*shift1,'r--', 'Linewidth', 2);  hold on
 %   c = plot(horizon, IRF_compl_w_2(K+i,1:h+1)*shift2,'b:', 'Linewidth', 2);  hold on
 %   plot(horizon,zeros(h+1,1), 'k')
 %   
 %   title_i = strcat(name(i)) ;
 %   title(title_i,'FontSize',font_size)
 %      ylabel(labelIRF(i),'FontSize',font_size)
 %   axis tight
 %   set(gca,'box','off','FontSize',font_size)
 %   
 %   if i == 4
 %       leg = legend([a b c], {'Baseline', 'Shock 1', 'Shock2'},'FontSize',9,'Location','northeast');
 %   end
%end

%% ----------------------------------------------------------------------- %
% Table 2 Output
% Hypothesis Tests
% Test 1
Wald1;
pval = mean(Wald1_r > Wald1);
% Test 2
Wald2;
pval_2 = mean(Wald2_r>Wald2);

%---------------------------------
% Table 3
% Inference on Lambda
Lambda_prct1 = prctile(Lambda_mat',[5 95],1);
Lambda_prct2 = prctile(Lambda_mat',[16 84],1);
Lambda_sign  = [Lambda_prct1; Lambda_prct2; diag(Lambda)'];
Lambda_sign([1:5],:) = Lambda_sign([1 3 5 4 2],:);

%---------------------------------
% Table 4
% test statistic for different Lambdas
wald_lmns = test_LMNS([T K],ann_use'+1,Uhat',[Sigma1_hat; Sigma0_hat],[eye(K); Lambda], p); 
test_lmns = wald_lmns.waldTestNK;  
%---------------------------------


%% ----------------------------------------
%  Run descriptive statistics
%-----------------------------------------
descriptive_results = descriptiveStats(wheatf, cornf, sunoil_f, soybean_f, rice_f, ...
    rapeseed_f, cocoa_f, coffee_f, sugar_f, index_own, spagri, oil_f, ...
    natural_gas_f,natgas_ttf_m, natgas_ttf_y, natgas_nbp_m, natgas_nbp_y, natgas_psv_m, natgas_psv_y, ...
    vstoxx, euronext, ger10ydiff, useuro, natgas_india, natgas_japan);
disp(descriptive_results);

%% ---------------------------------------
% ADF and PP test
% Specify lags to test
lags = [1, 2, 3, 4, 5];

% Perform ADF tests
adfPValueResults = performADFTests(lags, wheatf, cornf, sunoil_f, soybean_f, rice_f, ...
    rapeseed_f, cocoa_f, coffee_f, sugar_f, index_own, spagri, oil_f, ...
    natural_gas_f,natgas_ttf_m, ...
    vstoxx, euronext, ger10ydiff, useuro, natgas_india, natgas_japan, ...
    wheatfdif, cornfdif, sunoil_fdif, soybean_fdif, rice_fdif, ...
    rapeseed_fdif, cocoa_fdif, sugar_fdif,index_owndif, spagridif, ...
    oil_fdif, natural_gas_fdif, natgas_ttf_mdif, ...
    vstoxxdif, euronextdif, ger10ydif, useurodif, natgas_indiadif, natgas_japandif);

% Display the p-value results
disp('ADF Test P-Values:');
disp(adfPValueResults);

% Perform PP tests
ppPValueResults = performPPTests(lags, wheatf, cornf, sunoil_f, soybean_f, rice_f, ...
    rapeseed_f, cocoa_f, coffee_f, sugar_f, index_own, spagri, oil_f, ...
    natural_gas_f,natgas_ttf_m, ...
    vstoxx, euronext, ger10ydiff, useuro, natgas_india, natgas_japan, ...
    wheatfdif, cornfdif, sunoil_fdif, soybean_fdif, rice_fdif, ...
    rapeseed_fdif, cocoa_fdif, sugar_fdif,index_owndif, spagridif, ...
    oil_fdif, natural_gas_fdif, natgas_ttf_mdif, ...
    vstoxxdif, euronextdif, ger10ydif, useurodif, natgas_indiadif, natgas_japandif);

% Display the p-value results
disp('PP Test P-Values:');
disp(ppPValueResults);


%% Function that outputs descriptive statistics
function statsTable = descriptiveStats(varargin)
    % Input:
    % varargin: Multiple input data series
    % Output:
    % statsTable: Table containing mean, median, standard deviation, range, skewness, kurtosis for each data series

    % Check number of input arguments
    numDataSeries = nargin;

    if numDataSeries == 0
        error('At least one data series is required.');
    end

    % Preallocate variables for statistics
    means = zeros(numDataSeries, 1);
    medians = zeros(numDataSeries, 1);
    stdDevs = zeros(numDataSeries, 1);
    ranges = zeros(numDataSeries, 1);
    skews = zeros(numDataSeries, 1);
    kurts = zeros(numDataSeries, 1);
    names = cell(numDataSeries, 1);

    % Store varargin in a temporary variable
    dataSeriesCell = varargin;

    % Calculate statistics for each data series
    for i = 1:numDataSeries
        currentData = dataSeriesCell{i};
        names{i} = inputname(i); % Get the name of the variable passed to the function

        % Replace Inf and NaN values with the last valid observation
        currentData = fillmissing(currentData, 'previous');

        if isempty(names{i})
            names{i} = sprintf('DataSeries%d', i);
        end

        if ~isvector(currentData)
            error('Each data series must be a vector.');
        end

        means(i) = mean(currentData);
        medians(i) = median(currentData);
        stdDevs(i) = std(currentData);
        ranges(i) = range(currentData);
        skews(i) = skewness(currentData);
        kurts(i) = kurtosis(currentData);

    end

    % Create a table with all statistics
    statsTable = table(means, medians, stdDevs, ranges, skews, kurts, ...
        'VariableNames', {'Mean', 'Median', 'StdDev', 'Range', 'Skewness', 'Kurtosis'}, ...
        'RowNames', names);

end

function adfPValues = performADFTests(lags, varargin)
    % varargin now directly takes the time series as separate inputs

    numSeries = numel(varargin);
    numLags = numel(lags);

    % Initialize result table
    adfPValues = cell(numSeries, numLags);
    seriesNames = cell(numSeries, 1);

    % Get series names and perform ADF test
    for i = 1:numSeries
        seriesNames{i} = inputname(i + 1); % Offset by 1 due to lags argument
        currentData = varargin{i};
        for j = 1:numLags
            % Perform ADF test
            [~, pValue, ~, ~, ~] = adftest(currentData, 'lags', lags(j));
            % Store p-value
            adfPValues{i, j} = pValue;
        end
    end

    % Convert results to table for better readability
    adfPValues = cell2table(adfPValues, 'RowNames', seriesNames, ...
                                        'VariableNames', arrayfun(@(x) ['Lag' num2str(x)], lags, 'UniformOutput', false));
end

function ppPValues = performPPTests(lags, varargin)
    % varargin now directly takes the time series as separate inputs

    numSeries = numel(varargin);
    numLags = numel(lags);

    % Initialize result table
    ppPValues = cell(numSeries, numLags);
    seriesNames = cell(numSeries, 1);

    % Get series names and perform PP test
    for i = 1:numSeries
        seriesNames{i} = inputname(i + 1); % Offset by 1 due to lags argument
        currentData = varargin{i};
        for j = 1:numLags
            % Perform PP test
            [~, pValue, ~, ~, ~] = pptest(currentData, 'lags', lags(j));
            % Store p-value
            ppPValues{i, j} = pValue;
        end
    end

    % Convert results to table for better readability
    ppPValues = cell2table(ppPValues, 'RowNames', seriesNames, ...
                                        'VariableNames', arrayfun(@(x) ['Lag' num2str(x)], lags, 'UniformOutput', false));
end

