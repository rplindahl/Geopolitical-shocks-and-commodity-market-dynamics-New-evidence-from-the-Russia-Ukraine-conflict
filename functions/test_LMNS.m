 function obj = test_LMNS(T,regimeRealisations,residuals,Sigma,Lambda,lags)
 
 % test_LMNS([T K],ann_use'+1,Uhat',[Sigma1_hat; Sigma0_hat],[eye(K); Lambda], p);
 
 % from Lütkepohl, Meitz, Netsunajev and Saikkonen Econometrics Journal 2021
 % need GetSigmaForRegime.m
 % T: [T k]
 % regimeRealisations: e.g., [1 2 1 1 2]' vector that gives regimes for each t
 % residuals: TxK
 % Sigma: [Sigma1; Sigma2] (2K x K)
 % Lambda: [eye(K) Lambdas]: lambdas for first and second regime stacked
 % lags: VAR lags specified
 
            test_count = 1;
            for s = 0:T(1,2)-1
                for r = 2:T(1,2)-s
                    z_m = zeros(T(1,2), max(regimeRealisations));
                    w_m = zeros(T(1,2), max(regimeRealisations));
                    kappa = zeros ( max(regimeRealisations), 1);
                    for regime = 1 : max(regimeRealisations)
                        for k = 1: T(1,2)
                            z_m(k, regime) = sum((residuals(regimeRealisations(:, 1) == regime, k) - ...
                            mean(residuals(regimeRealisations(:, 1)== regime, k))) .^4);
                            sigma = GetSigmaForRegime(Sigma, T, regime); 
                            z_m(k, regime) = z_m(k, regime) - 6 * sigma(k,k)^2;
                            z_m(k, regime) = z_m(k, regime)/(sum(regimeRealisations(:, 1) == regime) - 4);

                            w_m_t = (sum(regimeRealisations(:, 1) == regime))/(sum(regimeRealisations(:, 1) == regime) - 1);
                            w_m(k, regime) = w_m_t*(sigma(k,k)^2 - z_m(k, regime)/sum(regimeRealisations(:, 1) == regime));
                        end
                        kappa(regime, 1) = 1/(3 * T(1,2)) * sum(z_m(:, regime) ./ w_m(:, regime) ) - 1;
                    end
                    tau = sum(regimeRealisations(:, 1) == 1) / (T(1,1) - lags);
                    c_tau2 = ( (1 + kappa(1,1)) / tau + ( (1 + kappa(2,1)) /(1-tau)))^-1;
                    c_tau2_nk = ( (1) / tau + ( (1) /(1-tau)))^-1;
                    lambda_sum = 0;
                    lambda2 = GetSigmaForRegime(Lambda, T, 2);
                    lambda_sum_r = 0;
                    for n = s+1 : s + r
                        lambda_sum = lambda_sum + log(lambda2(n, n));
                        lambda_sum_r = lambda_sum_r + lambda2(n,n);
                    end
                    lambda_sum_r = log(lambda_sum_r / r);
                    Qr = c_tau2 * (-1*(T(1,1)-lags) * lambda_sum + (T(1,1)-lags) * r * lambda_sum_r );
                    Qr_nk = c_tau2_nk * (-1*(T(1,1)-lags) * lambda_sum + (T(1,1)-lags) * r * lambda_sum_r );
                    degrees_of_freedom = 0.5 * (r+2) * (r-1);

                    Qr_star = 2 * Qr / ((r+2)*(r-1));
                    p_value_star = 1 - fcdf(Qr_star, 0.5*(r+2)*(r-1), (T(1,1)-lags - T(1,2)*lags - 1));
                    
                    p_value = 1 - chi2cdf(Qr, degrees_of_freedom);
                    p_value_nk = 1 - chi2cdf(Qr_nk, degrees_of_freedom);
                    
                    obj.waldTest(test_count, 1) = s;
                    obj.waldTest(test_count, 2) = r;
                    obj.waldTest(test_count, 3) = Qr;
                    obj.waldTest(test_count, 4) = degrees_of_freedom;
                    obj.waldTest(test_count, 5) = p_value;

                    obj.waldTest(test_count, 6) = Qr_star;
                    obj.waldTest(test_count, 7) = p_value_star;
                    
                    
                    obj.waldTestNK(test_count, 1) = s;
                    obj.waldTestNK(test_count, 2) = r;
                    obj.waldTestNK(test_count, 3) = Qr_nk;
                    obj.waldTestNK(test_count, 4) = degrees_of_freedom;
                    obj.waldTestNK(test_count, 5) = p_value_nk;                    
                    test_count = test_count + 1;
                end
            end
        end