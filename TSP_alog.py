'''
引用包
'''
import concurrent.futures

import numpy as np
import random
import math
import pandas as pd
import os
import natsort
import threading
from time import ctime,sleep

'''
算法函数
'''

#最小插入
def Insert(D,i):
    n = D.shape[0]
    M=1000000
    row, col = np.diag_indices_from(D)
    D[row,col] = M
    j = D.shape[0]-1
    Circle=[j]
    schedule = [i for i in range(n)]
    schedule.pop(j)
    mincol,index = D[j,schedule].min(),D[j,schedule].argmin()
    j=schedule[index]
    Circle.append(j)
    schedule.pop(index)
    while len(Circle)<n:
        row_idx = np.array(Circle)
        col_idx = np.array(schedule)
        A = D[row_idx[:, None], col_idx]
        xybound = np.where(A==np.min(A).min())
        x,y = xybound[0],xybound[1]
        j=schedule[y[0]]
        insertcost=M
        for k in range(len(Circle)-1):
            testcost=D[Circle[k],j]+D[j,Circle[k+1]]-D[Circle[k],Circle[k+1]]
            if testcost<insertcost:
                insertcost=testcost;
                insertlocation=k;
        Circle = np.concatenate((Circle[:insertlocation+1],[j],Circle[insertlocation+1:]))
        Circle = Circle.tolist()
        schedule.pop(y[0])
    return np.array(Circle)

#临近插入
def Adjacent(D,i):
    n = D.shape[0]
    M = 1000000
    row, col = np.diag_indices_from(D)
    D[row,col] = M
    j = D.shape[0]-1
    Circle=[j]
    schedule=[i for i in range(n)]
    schedule.pop(j)
    while len(Circle)<n:
        mincol,index = D[j,schedule].min(),D[j,schedule].argmin()
        j=schedule[index]
        Circle.append(j)
        schedule.pop(index)
    return np.array(Circle)

#初始解矩阵
def InitPop(D,N):
    Chrom = np.zeros((N,N))
    for i in range(N):
        rand = np.random.uniform(0,1)
        if rand<=0.5:
            Circle = Insert(D,i)
        else:
            Circle = Adjacent(D,i)
        Chrom[i,:] = Circle
    return Chrom

#路径长度
def PathLength(D,Chrom):
    row,col = D.shape[0],D.shape[1]
    NIND = len(Chrom)
    L =  np.zeros((NIND,1))
    if math.sqrt(Chrom.size)!=Chrom.shape[0]:
        Chrom = Chrom.reshape(1,-1)
        NIND = 1
    for i in range(NIND):
        N = D.shape[0]
        p = Chrom[i,:].astype(np.int64)
        i1 = p[0:-1]
        i2 = p[1:]
        L[i,0] = np.sum(np.array(D.T).flatten()[(i2*row+i1).astype(np.int64).tolist()])+D[p[-1],p[0]]
    return L

#适应度函数
def Fitness(length):
    FitnV = 1./length
    return FitnV

#轮盘
def Sus(FitnV,Nsel):
    Nind,ans = FitnV.shape[0],FitnV.shape[1]
    cumfit = np.cumsum(FitnV)
    trials = cumfit[Nind-1]/Nsel*(random.uniform(0,1)+np.array([i for i in range(int(Nsel))]).T)
    Mf = np.repeat(cumfit,Nsel).reshape(len(cumfit),Nsel)
    Mt = np.repeat(trials,Nind).reshape(len(trials),Nind).T
    bound = np.where((Mt<Mf)&(np.concatenate([np.zeros((1,Nsel)),Mf[:Nind-1,:]])<=Mt))
    NewChrIx,ans = bound[0],bound[1]
    randMatrix = np.random.uniform(0,1,size=(Nsel,1))
    ans,shuf = np.sort(randMatrix),np.argsort(randMatrix)
    NewChrIx=NewChrIx[shuf]
    return NewChrIx

#选择基因
def Select(Chrom,FitnV,GGAP):
    NIND = Chrom.shape[0]
    Nsel = int(max(np.floor(NIND*GGAP+.5),2))
    ChrIx = Sus(FitnV,Nsel)
    SelCh = Chrom[ChrIx,:].reshape(Nsel,NIND)
    return SelCh

#基因交叉
def intercross(a, b):
    L = len(a)
    r1 = np.random.randint(L)
    r2 = np.random.randint(L)
    if r1 != r2:
        a0 = a
        b0 = b
        s, e = min(r1, r2), max(r1, r2)
        for i in range(s, e):
            a1 = a
            b1 = b
            a[i] = b0[i]
            b[i] = a0[i]
            x, y = np.where(a == a[i])[0], np.where(b == b[i])[0]
            i1 = x[x != i]
            i2 = y[y != i]
            if i1.size != 0:
                a[i1] = a1[i]
            if i2.size != 0:
                b[i1] = b1[i]
    return a, b

#重组基因
def Recombin(SelCh,Pc):
    Nsel = SelCh.shape[0]
    for i in range(0,Nsel-np.mod(Nsel,2),2):
        rand = np.random.uniform(0,1)
        if Pc>= rand:
            SelCh[i,:],SelCh[i+1,:]=intercross(SelCh[i,:],SelCh[i+1,:])
    return SelCh


#基因变异
def Mutate(SelCh,Pm):
    NSel,L = SelCh.shape[0],SelCh.shape[1]
    for i in range(NSel):
        rand = np.random.uniform(0,1)
        if Pm>=rand:
            R = np.array([i for i in range(L)])
            np.random.shuffle(R)
            SelCh[i,R[:2]] = SelCh[i,R[1::-1]]
    return SelCh

#基因反转/加速变异用
def Reverse(SelCh,D):
    row,col = SelCh.shape[0],SelCh.shape[1]
    ObjV = PathLength(D,SelCh)
    SelCh1 = SelCh
    for i in range(row):
        r1 = np.random.randint(col)
        r2 = np.random.randint(col)
        mininverse = min(r1,r2)
        maxinverse = max(r1,r2)
        SelCh1[i,mininverse:maxinverse] = SelCh1[i,mininverse:maxinverse][::-1]
    ObjV1 = PathLength(D,SelCh1)
    index = ObjV1<ObjV
    SelCh[index.flatten(),:] = SelCh1[index.flatten(),:]
    return SelCh

#更新种群
def Reins(Chrom,SelCh,ObjV):
    NIND = Chrom.shape[0]
    NSel = SelCh.shape[0]
    TobjV,index = np.sort(ObjV),np.argsort(ObjV,axis=0)
    Chrom = np.array(np.concatenate([Chrom[index[:(NIND-NSel),:].flatten()],SelCh]).tolist())
    return Chrom


#输出最优解
def OutputPath(R):
    N = len(R)
    a = int(np.where(R==0)[0])
    RR = np.concatenate([R[a:],R[0:a]])
    for i in range(N):
        if i!=N-1:
            print("{0}->".format(RR[i]),end="")
        else:
            print("{0}".format(RR[i]),end="")


#输出文件
def tocsv(R,Idpath,pathlength):
    N = len(R)
    a = int(np.where(R == 0)[0])
    RR = np.concatenate([R[a:], R[0:a]])
    idFile = pd.read_csv(Idpath,header=None).to_numpy().astype(np.int64)
    idx = RR
    listrange = np.array([i for i in range(N)])
    pathLen = np.concatenate([pathlength,[0 for i in range(N-1)]])
    uid = np.concatenate([[0],np.array(idFile[np.array(idx[1:]).astype(np.int64) - 1]).flatten()])
    resultdict = {'idx':idx,'uid':uid,'range':listrange,'pathlength':pathLen}
    result = pd.DataFrame(data=resultdict)
    result.to_csv(Idpath.split('_userId.csv')[0]+'_result.csv',encoding='utf_8_sig',index=None)
    return result

def tolist(R,Id_Dataframe,pathlength):
    N = len(R)
    a = int(np.where(R == 0)[0])
    RR = np.concatenate([R[a:], R[0:a]])
    idFile = Id_Dataframe
    idx = RR
    uid = np.concatenate([[0], np.array(idFile[np.array(idx[1:]).astype(np.int64) - 1]).flatten()])
    return idx,uid,pathlength


#多线程测试
def threadextra(formatedPath,userIdPath):
    Formated_fileList = natsort.natsorted(os.listdir(formatedPath))
    Userid_fileList = natsort.natsorted(os.listdir(userIdPath))
    print(Formated_fileList)
    print(Userid_fileList)
    if len(Formated_fileList)!=len(Userid_fileList):
        print('false')
        return
    else:
        threads = []
        for i in range(len(Formated_fileList)):
            threads.append(threading.Thread(target=TSPmain,args=(formatedPath+'/'+Formated_fileList[i],userIdPath+'/'+Userid_fileList[i])))
    return threads

#多进程
def concurrent_mutiprocess(formatedPath,userIdPath,workers):
    Formated_fileList = natsort.natsorted(os.listdir(formatedPath))
    Userid_fileList = natsort.natsorted(os.listdir(userIdPath))
    print(Formated_fileList)
    print(Userid_fileList)
    if len(Formated_fileList) != len(Userid_fileList):
        print('false')
        return
    futures = []
    result = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for i in range(len(Formated_fileList)):
            futures.append(executor.submit(TSPmain, formatedPath+'/'+Formated_fileList[i],userIdPath+'/'+Userid_fileList[i]))
        for future in concurrent.futures.as_completed(futures):
            result.append(future.result())

    return result




#主函数/多线程/多进程
def TSPmain(Dpath,Idpath):
    #完整用户矩阵
    # D = pd.read_csv(
    #     r'C:\Users\chenz\Documents\Tencent Files\1294191169\FileRecv\ChaiFenData (2)\Distance_Mat\巴南2-距离矩阵.csv',
    #     header=None).iloc[1:,1:].to_numpy().astype(np.int64)
    #子矩阵
    D = pd.read_csv(Dpath,header=None).to_numpy().astype(np.int64)
    MAXGEN = 2000
    Pc = 0.8
    Pm = 0.15
    GGAP = 0.7
    N = D.shape[0]
    Chrom = InitPop(D, N)
    print('初始种群中的一个随机值:')
    OutputPath(Chrom[1, :])
    Rlength = PathLength(D, Chrom[1, :])
    print('')
    print('总距离:{0}'.format(Rlength[0][0]))
    gen = 0
    while gen < MAXGEN:
        ObjV = PathLength(D, Chrom)
        FitnV = Fitness(ObjV)
        SelCh = Select(Chrom, FitnV, GGAP)
        SelCh = Recombin(SelCh, Pc)
        SelCh = Mutate(SelCh, Pm)
        SelCh = Reverse(SelCh, D)
        Chrom = Reins(Chrom, SelCh, ObjV)
        gen = gen + 1
    ObjV = PathLength(D, Chrom)
    minObjV, minInd = np.min(ObjV), np.argmin(ObjV)
    result = tocsv(Chrom[minInd,:], Idpath.to_numpy().astype(np.int64),ObjV[minInd])
    return result

#单线程测试用
def Siglemain(Dpath,Idpath):
    # 完整用户矩阵
    # D = pd.read_csv(
    #     r'C:\Users\chenz\Documents\Tencent Files\1294191169\FileRecv\ChaiFenData (2)\Distance_Mat\巴南2-距离矩阵.csv',
    #     header=None).iloc[1:,1:].to_numpy().astype(np.int64)
    # 子矩阵
    D = pd.read_csv(Dpath,header=None).to_numpy().astype(np.int64)
    Idpath = Idpath
    MAXGEN = 2000
    Pc = 0.8
    Pm = 0.15
    GGAP = 0.7
    N = D.shape[0]
    Chrom = InitPop(D, N)
    print('初始种群中的一个随机值:')
    OutputPath(Chrom[1, :])
    Rlength = PathLength(D, Chrom[1, :])
    print('')
    print('总距离:{0}'.format(Rlength[0][0]))
    gen = 0
    while gen < MAXGEN:
        ObjV = PathLength(D, Chrom)
        FitnV = Fitness(ObjV)
        SelCh = Select(Chrom, FitnV, GGAP)
        SelCh = Recombin(SelCh, Pc)
        SelCh = Mutate(SelCh, Pm)
        SelCh = Reverse(SelCh, D)
        Chrom = Reins(Chrom, SelCh, ObjV)
        gen = gen + 1
    ObjV = PathLength(D, Chrom)
    minObjV, minInd = np.min(ObjV), np.argmin(ObjV)
    # 输出最优解的路线和总距离
    print('最优解：')
    tocsv(Chrom[minInd, :], Idpath,ObjV[minInd])
    OutputPath(Chrom[minInd, :])
    print('')
    print('总距离:{0}'.format(ObjV[minInd]))

#输入两个pandas，输出一个顺序列表，一个用户id列表，一个最优距离
def SiglemainDataFrame(D_DataFrame,Id_DataFrame):
    # 完整用户矩阵
    # D = pd.read_csv(
    #     r'C:\Users\chenz\Documents\Tencent Files\1294191169\FileRecv\ChaiFenData (2)\Distance_Mat\巴南2-距离矩阵.csv',
    #     header=None).iloc[1:,1:].to_numpy().astype(np.int64)
    # 子矩阵
    D = D_DataFrame.iloc[1:,1:].to_numpy().astype(np.int64)
    Id = Id_DataFrame.to_numpy().astype(np.int64)
    MAXGEN = 2000
    Pc = 0.8
    Pm = 0.15
    GGAP = 0.7
    N = D.shape[0]
    Chrom = InitPop(D, N)
    print('初始种群中的一个随机值:')
    OutputPath(Chrom[1, :])
    Rlength = PathLength(D, Chrom[1, :])
    print('')
    print('总距离:{0}'.format(Rlength[0][0]))
    gen = 0
    while gen < MAXGEN:
        ObjV = PathLength(D, Chrom)
        FitnV = Fitness(ObjV)
        SelCh = Select(Chrom, FitnV, GGAP)
        SelCh = Recombin(SelCh, Pc)
        SelCh = Mutate(SelCh, Pm)
        SelCh = Reverse(SelCh, D)
        Chrom = Reins(Chrom, SelCh, ObjV)
        gen = gen + 1
    ObjV = PathLength(D, Chrom)
    minObjV, minInd = np.min(ObjV), np.argmin(ObjV)
    # 输出最优解的路线和总距离
    print('最优解：')
    #tocsv(Chrom[minInd, :], Idpath,ObjV[minInd])
    idx,uid,length = tolist(Chrom[minInd, :],Id,ObjV[minInd])
    print(idx)
    print(uid)
    print(length)
    OutputPath(Chrom[minInd, :])
    print('')
    print('总距离:{0}'.format(ObjV[minInd]))
    return uid,length


