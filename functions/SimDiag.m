function [B,Lambda] = SimDiag(Sigma1,Sigma2, sort)
%%%% fuction to simultaneously diagonalize two covariance matrices
%%%% computes B and Lambda such that Sigma1=B*B' and Sigma2=B*Lambda*B'
%%%%%% algorithm from (Golub/van Loan, Alg 8.7.1)
%%%%%%
if nargin() < 3
    sort = 0;
end
G = chol(Sigma1);
G=G';
C = inv(G)*Sigma2*inv(G)';
[U,T] = schur(C);           % T-s are also eigenvalues of the matrix C
Lambda=diag(diag(T));
B=G*U;

%% sort the L-s in decreasing order
if sort == 1
    L = length(Lambda);
    for i = 1: L-1
        for j=i+1:L
            if Lambda(i,i) < Lambda(j,j)
                L_temp = Lambda(i,i);
                Lambda(i,i) = Lambda(j,j);
                Lambda(j,j) = L_temp;

                B_temp(:,1) =  B(:,i);
                B(:,i) = B(:,j);
                B(:,j) = B_temp(:, 1);
            end
        end
    end
end

