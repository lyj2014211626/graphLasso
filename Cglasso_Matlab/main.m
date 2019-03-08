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
%     saveas(gcf,strcat(strcat('graph/',strName(ind,1:2)),'经验逆协方差.png'));
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
    excelName={'Name','双顶径', '大脑中动脉Vs', '大脑中动脉Vd', '大脑中动脉RI', '大脑中动脉PI', '脐动脉腹内Vs', '脐动脉腹内Vd',...
        '脐动脉腹内RI', '脐动脉腹内PI', '心轴（度）', '胸廓横径', '心脏横径', '胸廓面积', '心脏面积', '左房横径', '右房横径', '左室横径',...
        '右室横径', '主动脉', '主动脉峡部', '降主动脉1', '肺动脉', '左肺动脉', '右肺动脉', '动脉导管直径中间段', '二尖瓣 E峰', '二尖瓣 A峰',...
        '三尖瓣 E峰', '三尖瓣 A峰', '主动脉瓣', '肺动脉瓣', '主动脉弓Vs', '主动脉弓Vd', '动脉导管Vs', '动脉导管Vd', '静脉导管S峰', '静脉导管D峰',...
        '静脉导管A峰', '肺静脉S峰', '肺静脉D峰', '肺静脉A峰', '卵圆孔直径', '心率'};
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
    saveas(gcf,strcat(strcat('graph/',strName(ind,1:2)),'协方差.png'));
end
