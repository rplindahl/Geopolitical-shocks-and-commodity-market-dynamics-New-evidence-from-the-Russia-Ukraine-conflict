% Calculate c_w explicitly and do not set to 1, rather normalize b1 element, fmincon
% Calculation of first column of B matrix via event-observation heteroskedasticity
% /(Ti)


function [B1_base,B,fval,V_inv,S,V0,V1, eig_values] = ident_het(ann_chosen,T,K,Uhat,z0,A_con,b_con,Aeq,beq,lb,ub,nonlcon,options_con,ann_start)

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

eig_values = eig(Sigma1_hat/(Sigma0_hat));

% Estimate V0  and V1 the var-cov matrices of vech(Sigmai_hat)

% via formula on p.499 in Kilian,Lütkepohl (2017) and divided by T0 and T1
% to get a Wald test statistic
[V0,V1,~] = cov_sigmaU(Uhat,Uhat0, Uhat1, K);
V_inv = inv(V0+V1);
   

% minimize the objective function to obtain R
S = vech(Sigma1_hat - Sigma0_hat);


obj     = @(z) het_obj(z,S,(V0+V1),K);
[z_het, fval, ~, ~] = fmincon(obj,z0, A_con,b_con,Aeq,beq,lb,ub,nonlcon,options_con);  % quasi-newton method, exitflag 1: success, 0: solver stopped, 2: change in x smaller than step tolerance, 5:predicted decrease in obj fct was less than the function tolerance    
B1_base = z_het;
                      

B = [B1_base, zeros(size(B1_base,1),size(B1_base,1)-1)];

end

