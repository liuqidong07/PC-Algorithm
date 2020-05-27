# -*- encoding: utf-8 -*-
'''
@File    :   MyPC.py
@Time    :   2020/03/20 21:30:19
@Author  :   Liu Qidong
@Version :   1.0
@Contact :   dong_liuqi@163.com
'''

# here put the import lib
import numpy as np
import pandas as pd
from itertools import combinations
import math
from scipy.stats import norm
import networkx as nx
from matplotlib import pyplot as plt
from scipy.linalg import pinv2

#构建无向的skeleton
def make_skeleton(V, corr_matrix, sample_num):
    #输入节点列表
    node_num = len(V)
    C = np.ones((node_num, node_num))    #先构造一个完全图
    #构造空的分离集
    S = []
    for i in range(node_num):
        S.append([])
        for j in range(node_num):
            S[i].append([])
    
    pairs = []
    #构造点对列表
    for i in range(node_num):
        for j in range(node_num - i):
            if(i != (node_num - j - 1)):    #不要把(i,i)加进去
                pairs.append((i, (node_num - j - 1)))
            else:
                #在这里顺便把对角线的点变为0
                C[i, i] = 0
    
    l = -1    #初始化l
    while 1:
        l = l + 1
        flag = True    #用来判断是否所有的邻接集的大小都小于l
        for (i, j) in pairs:

            adj_set = get_adjSet(i, C, node_num)    #获得i的邻接点集
            

            if(C[i][j] == 1) & (len(adj_set) >= l):    #如果两个点邻接，且邻接点集大小大于l
                flag =False    #一旦进入这个if，说明存在，所以不应该推出循环
                adj_set.remove(j)    #邻接点集去掉j

                combin_set = combinations(adj_set, l)    #获得K集合长度为l的全部情况
                for K in combin_set:
                    if indep_judge(i, j, list(K), corr_matrix, sample_num):
                        #如果i,j在K的条件下独立，则去掉i,j间的边，并把K集合加入到分离集中
                        C[i][j] = 0
                        C[j][i] = 0
                        #pairs.remove((i, j))    #已经独立的点对不需要再进行判断了
                        S[i][j] = list(K)
                        S[j][i] = list(K)
                        break    #如果证明了i和j独立，则没必要再判断了
                    else:
                        continue
            else:
                continue
        
        #判断是否该推出循环
        if flag:
            break

    return C, S

#判断i和j是否在K集合的条件下相互独立
def indep_judge(i, j, K, corr_matrix, sample_num):
    #输入与判断独立性的i,j以及分离集K
    #corr_matrix是相关系数矩阵
    indep = True

    #计算条件偏相关系数
    if len(K) == 0:
        #条件集为空，那么条件偏相关系数就是两个点的相关系数
        r = corr_matrix[i, j]
    else:
        corr = corr_matrix[np.ix_([i] + [j] + K, [i] + [j] + K)]
        #先通过对相关系数矩阵求逆得到偏相关系数矩阵
        partial_corr = np.linalg.pinv(corr)    #求广义逆矩阵
        #partial_corr = pinv2(corr)
        #然后取偏相关系数求条件偏相关系数
        #条件偏相关系数求法：-1*pr(x,y)/sqrt(pr(x,x)*pr(y,y))
        r = (-1 * partial_corr[0, 1]) / (math.sqrt(abs(partial_corr[0, 0] * partial_corr[1, 1])))
    
    #把r压缩到(-1, 1)之中
    r = min(0.99999, max(r, -0.99999))     #大于1的都归到1，小于-1都归到-1

    #进行fisher变换
    z = 0.5 * math.log1p((2 * r) / (1 - r))    #因为log1p是ln(1+x)，所以要减1
    #进行fisher变换后z服从正态分布，把其变成标准正态变量
    z_standard = z * math.sqrt(sample_num - len(K) - 3)

    #假设检验变量落在置信区间中
    alpha = 0.005    #置信度
    if 2 * (1 - norm.cdf(abs(z_standard))) >= alpha:
        indep = True
    else:
        indep = False
    return indep

#把无向的skeleton扩展成CPDAG
def extend_CPDAG(V, C, S):
    #输入无向的skeleton-C，分离集S
    #V是顶点集，C是skeleton，S是分离集
    G = C
    node_num = len(V)

    #先找到所有的三元组，因为之后的规则都是用在三元组上的。
    #先寻找所有相邻的二元组，再去找三元组
    pairs = []
    for i in range(node_num):
        for j in range(node_num):
            if(i != j):    #除去对角元素
                if(C[i][j] == 1):
                    pairs.append((i, j))
    
    triples = []
    for (i, j) in pairs:
        for k in range(node_num):
            if(C[j][k] == 1) & (k != i):
                triples.append([i, j, k])
    
    #构造PDAG。对于i-j-k，若j不在分离集(S[i][k])中，且i,k不相邻,则：i -> j <- k
    for [i, j, k] in triples:
        if (G[i][j] == 1) & (G[j][i] == 1) & (G[k][j] == 1) & (G[j][k] == 1) & (G[i][k] == 0) & (
            G[k][i] == 0):    #保证是i-j-k
            if j not in S[i][k]:
                G[j][i] = 0
                G[j][k] = 0

    #rule1:对于i -> j - k，如果i,k不相邻，则：i -> j -> k  （因为i<-j->k已经全部找出来了)
    for [i, j, k] in triples:
        if (G[i][j] == 1) & (G[j][i] == 0) & (G[k][j] == 1) & (G[j][k] == 1) & (G[i][k] == 0) & (
            G[k][i] == 0):
            G[k][j] = 0

    #rule2:对于i -> j -> k，且:i - k，则：i -> k （很显然是为了不能构成环）
    for [i, j, k] in triples:
        if (G[i][j] == 1) & (G[j][i] == 0) & (G[j][k] == 1) & (G[k][j] == 0) & (G[i][k] == 1) & (
            G[k][i] == 1):
            G[k][i] = 0

    #rule3:对于i-j->k,i-l->k,i-k，则：i->k
    for [i, j ,k] in triples:
        for [l, m, n] in triples:
            if (i == l) & (k == n):    #先找出i-j-k和i-l-k
                if (G[i][j] == 1) & (G[j][i] == 1) & (G[i][m] == 1) & (G[m][i] == 1) & (G[j][k] == 1
                ) & (G[k][j] == 0) & (G[m][k] == 1) & (G[k][m] == 0) & (G[i][k] == 1) & (G[k][i] == 1):
                    G[k][i] = 0

    #rule4:对于i-j->k,j->k->l,i-l，则：i->l
    for [i, j ,k] in triples:
        for [l, m, n] in triples:
            if (j == l) & (k == m):    #先找出i-j-k和j-k-l
                if (G[i][j] == 1) & (G[j][i] == 1) & (G[j][k] == 1) & (G[k][j] == 0) & (G[k][n] == 1
                ) & (G[n][k] == 0) & (G[i][n] == 1) & (G[n][i] == 1):
                    G[n][i] = 0

    return G

#获取点i在图G中的临界点集
def get_adjSet(i, G, node_num):
    adj = []
    for j in range(node_num):
        if G[i][j] == 1:
            adj.append(j)
    return adj

#使用networkx画图
def Draw(DAG, nodes):
    G = nx.DiGraph()
    nodes = nodes.tolist()
    G.add_nodes_from(nodes)    #加入所有的点
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if DAG[i][j] == 1:
                G.add_edge(nodes[i], nodes[j])
    nx.draw(G, with_labels=True)
    plt.show()


def PC(data_path):
    df = pd.read_csv(data_path, index_col=False)    #不读取行序号

    #计算相关系数矩阵
    corr = df.corr().values
    #样本的个数
    sample_num = df.values.shape[0]

    node = df.columns.values

    skeleton, separate = make_skeleton(node, corr, sample_num)
    CPDAG = extend_CPDAG(node, skeleton, separate)
    #print(CPDAG)
    #Draw(CPDAG, node)
    return CPDAG


if __name__ == '__main__':
    #data_path = r'./data/test.csv'
    data_path = r'./data/SimulateData.csv'
    PC(data_path)
    



