%%%% Generate data
clear;
%%% my data
strName=['01.txt';'02.txt';'03.txt';'04.txt'];
lamda=[0.2 0.02 0.075 0.53];
lamdaName=['lamda11';'lamda12';'lamda13';'lamda14'];
for ind=1:4
    A=importdata(strcat('weeks/',strName(ind,:)));
    [row,col]=size(A);
    sigma=cov(A);
    p=col;
    n=row;
%     imagesc(inv(sigma));
%     colormap(gray);
%     colorbar;
%     caxis([-1 1]);
%     title('Empirical Inverse Covariance Matrix')
%     saveas(gcf,strcat(strcat('graph/',strName(ind,1:2)),'������Э����.png'));
    %%%
    %p = 100; n =2*p;duq
    
    %% Sparse Sigma
    SigTrue = toeplitz([1,0.4,zeros(1,p-2)]);
    w = eig(SigTrue); delta = (p*min(w)-max(w))/(1-p); SigTrue = SigTrue + delta*eye(p);
    
    %% Dense Sigma
    % SigTrue = ones(p)+eye(p); % Dense
    Y = mvnrnd(zeros(p,1),SigTrue,n)';
    %S = Y*Y'/n; % sample covariance matrix;
    S=sigma;
    rho =lamda(ind); Rho = rho*ones(p)-rho*eye(p); % shrinkage parameters;
    
    SigInit = S; % Initial Value;
    %%%%%%%%%%% ECM algorithm
    %[Sig_ecm,C_ecm,loglik_ecm] = CglassoECM(S,Rho,SigInit,1e-3,1e4);
    %%%%%%%%%%% Coordinate descent algorithm
    excelName={'Name','˫����', '�����ж���Vs', '�����ж���Vd', '�����ж���RI', '�����ж���PI', '�궯������Vs', '�궯������Vd',...
        '�궯������RI', '�궯������PI', '���ᣨ�ȣ�', '�����ᾶ', '����ᾶ', '�������', '�������', '�󷿺ᾶ', '�ҷ��ᾶ', '���Һᾶ',...
        '���Һᾶ', '������', '������Ͽ��', '��������1', '�ζ���', '��ζ���', '�ҷζ���', '��������ֱ���м��', '����� E��', '����� A��',...
        '����� E��', '����� A��', '��������', '�ζ�����', '��������Vs', '��������Vd', '��������Vs', '��������Vd', '��������S��', '��������D��',...
        '��������A��', '�ξ���S��', '�ξ���D��', '�ξ���A��', '��Բ��ֱ��', '����'};
    [Sig_cd,C_cd,loglik_cd] = CglassoCD(S,Rho,SigInit,1e-3,200,1e4);
    dlmwrite(strcat(strcat('precision_Matrix/',strName(ind,1:2)),'Precision.txt'),Sig_cd,'delimiter','\t','precision','%6.4f','newline','pc')
    xlswrite(strcat(strcat('precision_Matrix/',strName(ind,1:2)),'Precision.xls'),excelName,'Sheet1')
    xlswrite(strcat(strcat('precision_Matrix/',strName(ind,1:2)),'Precision.xls'),excelName','Sheet1')
    xlswrite(strcat(strcat('precision_Matrix/',strName(ind,1:2)),'Precision.xls'),Sig_cd,'Sheet1','B2')
    imagesc(C_cd);
    colormap(gray);
    colorbar;
    caxis([-1 1]);
    title('Covariance Matrix')
    %title('Preision Matrix')
    saveas(gcf,strcat(strcat('graph/',strName(ind,1:2)),'Э����.png'));
end
