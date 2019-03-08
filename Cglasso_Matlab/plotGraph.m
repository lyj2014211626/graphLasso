fileName={'01Precision.xls';'02Precision.xls';'03Precision.xls';'04Precision.xls'};
nodeName={'双顶径', '大脑中动脉Vs', '大脑中动脉Vd', '大脑中动脉RI', '大脑中动脉PI', '脐动脉腹内Vs', '脐动脉腹内Vd',...
    '脐动脉腹内RI', '脐动脉腹内PI', '心轴（度）', '胸廓横径', '心脏横径', '胸廓面积', '心脏面积', '左房横径', '右房横径', '左室横径',...
    '右室横径', '主动脉', '主动脉峡部', '降主动脉1', '肺动脉', '左肺动脉', '右肺动脉', '动脉导管直径中间段', '二尖瓣 E峰', '二尖瓣 A峰',...
    '三尖瓣 E峰', '三尖瓣 A峰', '主动脉瓣', '肺动脉瓣', '主动脉弓Vs', '主动脉弓Vd', '动脉导管Vs', '动脉导管Vd', '静脉导管S峰', '静脉导管D峰',...
    '静脉导管A峰', '肺静脉S峰', '肺静脉D峰', '肺静脉A峰', '卵圆孔直径', '心率'};
nodeNum={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26'...
    ,'27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43'};
nodenum=1:43;
diff1={'5','14','20','23','40'};
diff2={'1','4','11','15','17','19','43'};
for index=4:4
    adjMatrix=xlsread(strcat('precision_Matrix/',char(fileName(index,:))));%图的邻接矩阵
    adjMatrix(abs(adjMatrix)<1e-3)=0;%将小于一定值的置零
    %adjMatrix=adjMatrix-eye(size(adjMatrix,1));
    G=graph(adjMatrix,nodeNum,'upper','OmitSelfLoops');%'OmitSelfLoops'去掉自环的边 %'upper'用邻接矩阵的上三角
    edges=G.Edges;%图的边
    nodes=G.Nodes;
    emptyNode=find(sum(adjMatrix)-diag(adjMatrix)'==0);
    nodes(emptyNode,:)=[];
    G=rmnode(G,emptyNode);%删除度为0的节点
    
    deg = degree(G);%节点的度
    g=plot(G,'Layout','circle','LineWidth',1.8,'MarkerSize',10);%线宽和圆圈的大小
    for i=1:size(diff2,2)
        neigh=G.neighbors(diff2(i));%邻接节点
        highlight(g,diff2(i),'NodeColor','g')
        highlight(g,diff2(i),neigh,'EdgeColor','g')
    end
end