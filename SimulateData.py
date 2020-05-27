# -*- encoding: utf-8 -*-
'''
@File    :   SimulateData.py
@Time    :   2020/03/27 11:09:40
@Author  :   Liu Qidong
@Version :   1.0
@Contact :   dong_liuqi@163.com
'''

# here put the import lib

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

#np.random.seed(99)


def Simulate_Data(sparseness, varibles, data_num):
    '''
    模拟因果推断实验数据
    :param sparseness: 稀疏度，也是构造矩阵A时伯努利实验的概率
    :param varibles: 变量名称列表
    :param data_num: 欲构造的数据量
    :return:
    '''
    data = []
    varible_num = len(varibles)
    A = np.zeros((varible_num, varible_num))    #论文中所说的矩阵A
    B = np.zeros((varible_num, varible_num))    #矩阵B作为记录有向边的邻接矩阵
    #bernoulli_list = np.random.binomial(1, sparseness, ((varible_num - 1) * varible_num) // 2)    #第三个参数是求除了A矩阵下三角元素个数

    #把A的下三角进行伯努利试验赋值
    for i in range(1, varible_num):
        for j in range(i):
            B[j][i] = np.random.binomial(1, sparseness)
            if B[j][i] == 1:
                A[i][j] = np.random.uniform(0.9, 1)
                #B[i][j] = 1

    #利用矩阵A来构造数据
    for order in range(data_num):
        #normal_list = np.random.normal(0, 1, varible_num)    #构造服从正态分布的残差值
        X = np.zeros(varible_num)    #先初始化一组变量（结点）的值，方便后面使用
        X[0] = np.random.normal(0, 1)
        for i in range(1, varible_num):
            for k in range(i):
                X[i] += A[i][k] * X[k]
            X[i] += np.random.normal(0, 1)
        data.append(X)
    return np.array(data), B



def Data2CSV(labels, data, save_path):
    '''
    把模拟数据写文件
    :param labels: 结点名称的列表
    :param data: 数据
    :return: 无
    '''
    #save_path = r'./SimulateData.csv'
    pd_dict = {}
    for label in labels:
        pd_dict[label] = []
    #要把一条条的数据展开成一列一列的，这样才能构建DataFrame
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            pd_dict[labels[j]].append(data[i][j])
    df = pd.DataFrame(pd_dict)
    df.to_csv(save_path, index = None)



def Draw(matrix):
    '''
    使用networkx把模拟的DAG画出来
    :param matrix: 模拟的邻接矩阵A
    :return:
    '''
    G = nx.DiGraph()
    labels = [i for i in range(len(matrix))]
    G.add_nodes_from(labels)
    edges = []
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if(matrix[i][j] == 1):
                edges.append((i, j))
    G.add_edges_from(edges)
    nx.draw(G, with_labels=True)
    plt.savefig('./result/simulate.png')
    plt.show()


def simulate_data(n, p, EN, path):
    #n表示样本的数量，p是节点个数，EN=s(p-1)是每个节点的邻接节点个数的期望，s为稀疏度
    labels = []
    for i in range(p):
        labels.append(i)    #构建标签
    s = EN / (p - 1)
    data, adjacency_matrix = Simulate_Data(s, labels, n)
    Data2CSV(labels, data, path)
    #Draw(adjacency_matrix)
    return adjacency_matrix



if __name__ == '__main__':
    n = 1000
    p = 7
    EN = 2
    path = r'./SimulateData.csv'
    simulate_data(n, p, EN, path)
