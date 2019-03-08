fileName={'01Precision.xls';'02Precision.xls';'03Precision.xls';'04Precision.xls'};
nodeName={'˫����', '�����ж���Vs', '�����ж���Vd', '�����ж���RI', '�����ж���PI', '�궯������Vs', '�궯������Vd',...
    '�궯������RI', '�궯������PI', '���ᣨ�ȣ�', '�����ᾶ', '����ᾶ', '�������', '�������', '�󷿺ᾶ', '�ҷ��ᾶ', '���Һᾶ',...
    '���Һᾶ', '������', '������Ͽ��', '��������1', '�ζ���', '��ζ���', '�ҷζ���', '��������ֱ���м��', '����� E��', '����� A��',...
    '����� E��', '����� A��', '��������', '�ζ�����', '��������Vs', '��������Vd', '��������Vs', '��������Vd', '��������S��', '��������D��',...
    '��������A��', '�ξ���S��', '�ξ���D��', '�ξ���A��', '��Բ��ֱ��', '����'};
nodeNum={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26'...
    ,'27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43'};
nodenum=1:43;
diff1={'5','14','20','23','40'};
diff2={'1','4','11','15','17','19','43'};
for index=4:4
    adjMatrix=xlsread(strcat('precision_Matrix/',char(fileName(index,:))));%ͼ���ڽӾ���
    adjMatrix(abs(adjMatrix)<1e-3)=0;%��С��һ��ֵ������
    %adjMatrix=adjMatrix-eye(size(adjMatrix,1));
    G=graph(adjMatrix,nodeNum,'upper','OmitSelfLoops');%'OmitSelfLoops'ȥ���Ի��ı� %'upper'���ڽӾ����������
    edges=G.Edges;%ͼ�ı�
    nodes=G.Nodes;
    emptyNode=find(sum(adjMatrix)-diag(adjMatrix)'==0);
    nodes(emptyNode,:)=[];
    G=rmnode(G,emptyNode);%ɾ����Ϊ0�Ľڵ�
    
    deg = degree(G);%�ڵ�Ķ�
    g=plot(G,'Layout','circle','LineWidth',1.8,'MarkerSize',10);%�߿��ԲȦ�Ĵ�С
    for i=1:size(diff2,2)
        neigh=G.neighbors(diff2(i));%�ڽӽڵ�
        highlight(g,diff2(i),'NodeColor','g')
        highlight(g,diff2(i),neigh,'EdgeColor','g')
    end
end