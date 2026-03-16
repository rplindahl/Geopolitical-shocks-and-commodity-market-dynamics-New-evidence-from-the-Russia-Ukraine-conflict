% Test statistic H0: Sigma0 = Sigma 1 by using a random bootstrap
% assignment of announcement dates

function [F_alt] = test_sigma(Uhat,announcement_use_m, announcement_count, K)

% Resid (KxT)
% Estimate the var-cov matrices of  follows 


Uhat_0 = Uhat(announcement_use_m == 0);  % split the residuals, column vector now
Uhat_0 = reshape(Uhat_0,size(Uhat,1),size(Uhat,2)-announcement_count);  % shape into matrix of form like Resid with appropriate T
Uhat_1 = Uhat(announcement_use_m == 1);
Uhat_1 = reshape(Uhat_1, size(Uhat,1), announcement_count);

T0 = length(Uhat_0);                    
T1 = length(Uhat_1);

mean_a = 1/(T0+T1)*(Uhat*Uhat');


vech_u_sq_dif = zeros(K*(K+1)/2,K*(K+1)/2);
for i = 1:length(Uhat)
    sq_dif   = (vech( Uhat(:,i)*Uhat(:,i)' - mean_a )) * (vech( Uhat(:,i)*Uhat(:,i)' - mean_a ))';
    vech_u_sq_dif = sq_dif + vech_u_sq_dif;
end
W = vech_u_sq_dif/(T0+T1);


V     = W/T0 + W/T1;            % assume same vcov for the two sets in this test
% V_inv = inv(V);


Sigma1_hat = (Uhat_1*Uhat_1')/T1;
Sigma0_hat = (Uhat_0*Uhat_0')/T0;
S = vech(Sigma1_hat - Sigma0_hat);


F_alt = S' /V * S;  

    
end