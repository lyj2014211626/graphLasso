function [Sig,C,obj] = CglassoCD(Q,Rho,Sig,tol,outermax,innermax)
% CglassoCD solves the covariance graphical lasso problem (Bien and Tibshirani 2011, Biometrika):
% Sigma_cd  = argmin { logdet(Sigma)+trace(Q/Sig)+sum(sum(Rho.*abs(Sig))) }
% using coordinate descent algorithm (Wang 2012)
% % Input:
%   Q:  p*p sample covariance matrix  Q = Y*Y'/n;
%   Rho: p*p shrinkage parameter matrix
%   Sig: inital value
%   tol: threshold for the change of objective functions for terminating iterations
%   outermax: maximum number of itermations of the outer loop
%   innermax: maximum number of iterations of the inner loop for the lasso problem
%
%   Output:
%   Sig: clustering point of Sig
%   C:   inv(Sig)
%   obj: minimum value of the objective function:
%                 obj = logdet(Sigma)+trace(Q/Sig)+sum(sum(Rho.*abs(Sig)));


% Written by Hao Wang @ U of South Carolina

p = size(Q,1);%the number of rows about Q

%生成ind_noi_all用来选择生成块矩阵
ind_noi_all = zeros(p-1,p);
for i = 1:p
    ind_noi = [1:p];
    ind_noi(i) = [];
    ind_noi_all(:,i) = ind_noi;
end

%% Initial value
C = inv(Sig);
nedge_est = (sum(sum(abs(Sig)>1e-3))-p)/2;% the number of edges
obj = logdet(Sig)+trace(C*Q)+sum(sum(Rho.*abs(Sig)));%objective function
disp(['iter = 0 n_edge = ',num2str(nedge_est), ' obj = ',num2str(obj)]);

%%
optTol = 1e-3;%bias 误差
%外层循环 最大迭代次数
for iter1 = 1:outermax 
    obj_old = obj;%开始的目标值
    %循环次数为特征变量的个数
    for i = 1:p           
        ind_noi = ind_noi_all(:,i);       
        C11 = C(ind_noi,ind_noi); C12 = C(ind_noi,i);
        S11 = Q(ind_noi,ind_noi); S12 = Q(ind_noi,i);
        invSig11 = C11 - C12*C12'/C(i,i);  
        invSig11S12 = invSig11*S12;
        W1 = invSig11*S11*invSig11;
        
        %% Update gamma
        beta = Sig(ind_noi,i);
        b = Rho(i,i)/2; c = (beta'*W1*beta-2*beta'*invSig11S12+Q(i,i))/2;
        a = 1/2;
        if b == 0
            gam = c/a;
        else
            gam = (-a+sqrt(a^2+4*b*c))/(2*b);
        end
        
        %% Update beta
        V = W1/gam+Rho(i,i)*invSig11;
        u = invSig11S12./gam;
                
        %% below Pathwise coordinate descent
        Vbeta = V*beta;        
        for iter2 = 1:innermax           
            beta_old = beta;    
            for j = 1:p-1
                betaj_old = beta(j);
                x = u(j)-Vbeta(j)+V(j,j)*beta(j);%减去k==j的情况
                beta(j) = max(0,abs(x)-Rho(i,ind_noi(j)))*sign(x)/V(j,j);
                beta_change = beta(j)-betaj_old;
                if beta_change~=0
                    Vbeta = Vbeta+beta_change*V(:,j);
                end           
            end
            
            if max(abs(beta_old-beta))< optTol
                break;
            end           
        end
     
        Sig(i,ind_noi) = beta';
        Sig(ind_noi,i) = beta;
        Sig(i,i) = gam+beta'*invSig11*beta;

        %% Below updating Precision matrix according to one-column change of precision matrix
        invSig11beta = invSig11*beta;
        C(ind_noi,ind_noi) = invSig11+invSig11beta*invSig11beta'./gam;
        C12 = -invSig11beta/gam;
        C(ind_noi,i) = C12;
        C(i,ind_noi) = C12';
        C(i,i) = 1/gam;
               
    end
    
    obj = logdet(Sig)+trace(C*Q)+sum(sum(Rho.*abs(Sig)));  
    nedge_est = (sum(sum(abs(Sig)>1e-3))-p)/2;
    disp(['iter = ',int2str(iter1),' n_edge = ',num2str(nedge_est), ' obj = ',num2str(obj)])
      
    if (abs(obj-obj_old)<tol)
        break;
    end
    
end
