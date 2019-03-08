function [Sig,C,obj] = CglassoECM(Q,Rho,Sig,tol,itermax)
% CglassoECM solves the covariance graphical lasso problem (Bien and Tibshirani 2011, Biometrika):
% Sigma_cd  = argmin { logdet(Sigma)+trace(Q/Sig)+sum(sum(Rho.*abs(Sig))); }
% using ECM algorithm (Wang 2012)
% % Input:
%   Q:  p*p sample covariance matrix  Q = Y*Y'/n; 
%   Rho: p*p shrinkage parameter matrix
%   Sig: inital value
%   tol: threshold for the change of objective functions for terminating iterations
%   itermax: maximum number of iterations 
%   
%   Output:
%   Sig: clustering point of Sig
%   C:   inv(Sig)
%   obj: minimum value of the objective function: 
%              obj: logdet(Sigma)+trace(Q/Sig)+sum(sum(Rho.*abs(Sig)));


% Written by Hao Wang @ U of South Carolina

n = 1;
S = Q;
[p] = size(Q,1);

C = inv(Sig);

% indmx = reshape([1:p^2],p,p); 
% upperind = indmx(triu(indmx,1)>0); 
% 
% indmx_t = indmx';
% lowerind = indmx_t(triu(indmx_t,1)>0); 


Lambda = Rho.*n;

ind_noi_all = zeros(p-1,p);
for i = 1:p
       if i==1  
       ind_noi = [2:p]'; 
      elseif i==p
       ind_noi = [1:p-1]'; 
      else
       ind_noi = [1:i-1,i+1:p]';
       end
       
       ind_noi_all(:,i) = ind_noi;
end


% Initial Objective function
    obj = -logdet(C)+trace(C*Q)+sum(sum(Rho.*abs(Sig)));
    nedge_est=(sum(sum(abs(Sig)>1e-3))-p)/2;

    disp(['iter = 0 n_edge = ',num2str(nedge_est), ' obj = ',num2str(obj)]);

    

    for iter = 1:itermax
        
        
    % Current Sig and objective function    
       obj_old = obj;
       Sig_old = abs(Sig);
    

     

    for i = 1:p
        
      ind_noi = ind_noi_all(:,i);
      
       
      C11 = C(ind_noi,ind_noi); C12 = C(ind_noi,i);
      S11 = S(ind_noi,ind_noi); S12 = S(ind_noi,i); 
      
      Sig11 = Sig(ind_noi,ind_noi);
      invSig11 = C11 - C12*C12'/C(i,i);  invSig11S12 = invSig11*S12;    
      
      S11invSig11 = S11*invSig11;
      W1 = invSig11*S11invSig11;
      %% Update gamma
      %      pdf = gam^(-a)*exp(-b gam-c/gam);
 
      beta = Sig(ind_noi,i); 
      b = Lambda(i,i)/2; c = (beta'*W1*beta - 2*beta'*invSig11S12+S(i,i))/2;    
      a = n/2; 
      if b == 0
        gam = c/a;
      else
        gam = (-a+sqrt(a^2+4*b*c))/(2*b);  
      end
      
  
        
      %% Sample beta       
      A = (S11invSig11./gam+Lambda(i,i)*eye(p-1))\Sig11;
      B = diag(Sig_old(ind_noi,i)./Lambda(ind_noi,i));
      
      M = B/(A+B)*A;
      
      mu_i =  M*invSig11S12./gam;
      beta = mu_i;
        
              
        Sig(ind_noi,i) = beta;
        Sig(i,ind_noi) = beta;       
        Sig(i,i) = gam+beta'*invSig11*beta;
        
        
        %% Below updating Precision matrix according to one-column change of precision matrix
        invSig11beta = invSig11*beta;
        
        C(ind_noi,ind_noi) = invSig11+invSig11beta*invSig11beta'./gam;
        C12 = -invSig11beta/gam;
        C(ind_noi,i) = C12;
        C(i,ind_noi) = C12';
        C(i,i) = 1/gam;
        
            
    end
    
         obj =  -logdet(C)+trace(C*Q)+sum(sum(Rho.*abs(Sig)));

    nedge_est = (sum(sum(abs(Sig)>1e-3))-p)/2;
      disp(['iter = ',int2str(iter),' n_edge = ',num2str(nedge_est), ' obj = ',num2str(obj)])
         
         if( abs(obj-obj_old)<tol)
             break;
         end
         
    end