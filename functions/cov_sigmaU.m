% Calculation of the cov matrices of the variances vech(Sigmai_hat)
% denoted V0 and V1 via the formula from Lütkepohl & Kilian (2017) p.499


function [V0, V1, V_w] = cov_sigmaU(Uhat, Uhat_0, Uhat_1, K)

% Resid (KxT)
% Estimate the var-cov matrices of  follows 


T0 = length(Uhat_0);                    
T1 = length(Uhat_1);

mean_0 = 1/T0*(Uhat_0*Uhat_0');  % same as Sigma0_hat
mean_1 = 1/T1*(Uhat_1*Uhat_1');  % same as Sigma1_hat

vech_u_sq0_dif = zeros(K*(K+1)/2,K*(K+1)/2);
for i = 1:size(Uhat_0,2)
    sq_dif   = (vech( Uhat_0(:,i)*Uhat_0(:,i)' - mean_0 )) * (vech( Uhat_0(:,i)*Uhat_0(:,i)' - mean_0 ))';
    vech_u_sq0_dif = sq_dif + vech_u_sq0_dif;
end
W0 = vech_u_sq0_dif/T0;

vech_u_sq1_dif = zeros(K*(K+1)/2,K*(K+1)/2);
for i = 1:size(Uhat_1,2)
    sq_dif   = (vech( Uhat_1(:,i)*Uhat_1(:,i)' - mean_1 )) * (vech( Uhat_1(:,i)*Uhat_1(:,i)' - mean_1 ))';
    vech_u_sq1_dif = sq_dif + vech_u_sq1_dif ;
end
W1 = vech_u_sq1_dif/T1;

% V     = W0/T0 + W1/T1;
% V_inv = inv(V);


V0 =  W0/T0;            % Output 1
V1 = W1/T1;             % Output 2   



% Wald Test 1 Weighting matrix input (S=0) as in Wright (2012)
mean_all = 1/(T0+T1)*(Uhat*Uhat');  


vech_u_sq_dif = zeros(K*(K+1)/2,K*(K+1)/2);
for i = 1:size(Uhat,2)
    sq_dif   = (vech( Uhat(:,i)*Uhat(:,i)' - mean_all )) * (vech( Uhat(:,i)*Uhat(:,i)' - mean_all ))';
    vech_u_sq_dif = sq_dif + vech_u_sq_dif;
end
M = vech_u_sq_dif/(T0+T1);

V_w = (M/T0+M/T1);   % Wald Test 1 as in Wright (2012) Code


    
end