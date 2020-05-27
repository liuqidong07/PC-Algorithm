# -*- encoding: utf-8 -*-
'''
@File    :   evaluation.py
@Time    :   2020/03/28 18:32:32
@Author  :   Liu Qidong
@Version :   1.0
@Contact :   dong_liuqi@163.com
'''

# here put the import lib
from MyPC import PC
from SimulateData import simulate_data


def evaluation():
    n = 1000    #样本数量
    p = 40    #节点个数
    EN = 5    #邻接节点平均个数（期望）
    replications = 40

    path = r'./data/SimulateData.csv'

    SHD_mean = 0
    #进行重复实验取平均值
    for l in range(replications):
        true_adj = simulate_data(n, p, EN, path)
        predict_adj = PC(path)
        SHD = 0    #Structural Hamming Distance
        for i in range(len(true_adj)):
            for j in range(len(true_adj[0])):
                if true_adj[i][j] != predict_adj[i][j]:
                    SHD = SHD + 1
        SHD_mean += SHD/replications
        print('The %d replication SHD is %d' % (l, SHD))

    print('Structural Hamming Distance is:%d' % SHD_mean)


if __name__ == '__main__':
    evaluation()

