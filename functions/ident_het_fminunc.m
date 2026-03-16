% Fminunc to estimate b1, c_w set to 1
% Calculation of first column of B matrix via event-observation heteroskedasticity
% b1 gives position of event that first Residual belongs to


function [B1,B,B1_base, fval,V_inv,S,V0,V1] = ident_het_fminunc(ann_chosen,T,K,Uhat,x0,options,normalizer,pos_norm,ann_start)

% Create appropriate index variable for announcement days
% choose announcement type here
ann_use    = ann_chosen(ann_start:T,:)';             

% Create the two distinct var-cov matrices
Uhat_t  = Uhat';
Uhat1_t = Uhat_t(ann_use==1,:);  
Uhat0_t = Uhat_t(ann_use==0,:);
Uhat0   = Uhat0_t';
Uhat1   = Uhat1_t';

T0 = length(Uhat0);                    
T1 = length(Uhat1);

Sigma1_hat = (Uhat1*Uhat1')/(T1);
Sigma0_hat = (Uhat0*Uhat0')/(T0);


% Estimate V0  and V1 the var-cov matrices of vech(Sigmai_hat)

% via formula on p.499 in Kilian,Lütkepohl (2017) and divided by T0 and T1
% to get a Wald test statistic
[V0,V1,~] = cov_sigmaU(Uhat,Uhat0, Uhat1, K);
V_inv = inv(V0+V1);
   

% minimize the objective function to obtain R
S = vech(Sigma1_hat - Sigma0_hat);

% 1a: 6 vars, use V_inv, fminunc
fun     = @(x)(S - vech(x*x'))'/(V0+V1)*(S - vech(x*x'));
[x, fval, ~, ~] = fminunc(fun,x0, options);  % quasi-newton method, exitflag 1: success, 0: solver stopped, 2: change in x smaller than step tolerance, 5:predicted decrease in obj fct was less than the function tolerance    
B1_base = x;
                      

%% Compute structural impulse responses
% x = R1 now, u_t = Sum_(1 to k) (Ri * w_t), w_t structural errors
% we only have R1


norm = x(pos_norm)/normalizer;          
x_norm = x/norm;                           

B1 = x_norm;
% save('R1.mat','R1'); 
B = [B1, zeros(size(B1,1),size(B1,1)-1)];

end

