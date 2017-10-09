
# coding: utf-8

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys


# In[45]:

def group_data(data,tol=5):
    # group consecutive integers (data) and tolerate gaps of (tol)
    out = []
    last = data[0]
    for x in data:
        if x-last > tol:
            yield out
            out = []
        out.append(x)
        last = x
    yield out


def getAPDx(AP,t,percentile=0.95,x=0.9,relax=100,tol=0.5,check2groups=True,groupTol=5):
    ## return AP duration at x of AP amplitude, APDx, along with MDP, Peak, and APA
    numIdx = int(percentile*len(AP))
    # maximum voltage with given percentile
    maxSig = np.sort(AP)[numIdx]
    # AP amplitude
    APA = maxSig-AP[0]
    # voltage at x of APA
    relativeMaxSigx = ((1.0-x)*APA + AP[0])
    # first #relax closest points from relativeMaxSigx
    argmin_relax = np.argsort(np.abs(AP-relativeMaxSigx))[:relax]
    # check only two groups of such closest points (else cannot def APD)
    if check2groups:
        temp_argmin_relax = list(argmin_relax)
        temp_argmin_relax.sort()
        if len(list(group_data(temp_argmin_relax,groupTol)))>2:
            return np.inf,np.inf,np.inf,np.inf
    isDuration = np.abs(argmin_relax-argmin_relax[0])>relax
    if np.any(isDuration):
        getDurIdx = np.arange(relax)[isDuration][0]
        if np.abs(AP[argmin_relax[getDurIdx]]-relativeMaxSigx)<tol:
            APDx = np.abs(t[argmin_relax[0]] - t[argmin_relax[getDurIdx]])
            # maximum diastolic potential (MDP)
            MDP = np.min(AP)
            return APDx, MDP, maxSig, APA, relativeMaxSigx
        else:
            return np.inf,np.inf,np.inf,np.inf,np.inf
    else:
        return np.inf,np.inf,np.inf,np.inf,np.inf


def getAPar(AP,t,APDx,maxSig,relativeMaxSigx):
    ## return AP area ratio defined in Wang et al 2015
    recArea = (maxSig-relativeMaxSigx)*APDx
    APArea = 0
    i = 0
    while AP[i] < relativeMaxSigx:
        i += 1
    while AP[i] > relativeMaxSigx and AP[i] < maxSig:
        APArea += (AP[i]-relativeMaxSigx)*(t[i]-t[i-1])
        i += 1
    while AP[i] > maxSig:
        APArea += (maxSig-relativeMaxSigx)*(t[i]-t[i-1])
        i += 1
    while AP[i] > relativeMaxSigx:
        APArea += (AP[i]-relativeMaxSigx)*(t[i]-t[i-1])
        i += 1
    return float(APArea)/float(recArea)


def getSmoothAP(AP,t,t0=None,t1=None,t2=None):
    # time of t0:upstroke; t1=maxsignal; t2:baseline
    if t0==None or t1==None or t2==None:
        #do something else
        pass
    idx0 = np.where(t==t0)[0][0]
    idx1 = np.where(t==t1)[0][0]
    idx2 = np.where(t==t2)[0][0]
    lin = AP[idx0:idx1]
    qua = AP[idx1:idx2]
    fitlin = np.polyfit(t[idx0:idx1], lin, 1)
    fitqua = np.polyfit(t[idx1:idx2], qua, 4)
    plin = np.poly1d(fitlin)
    pqua = np.poly1d(fitqua)
    APlin = plin(t[idx0:idx1])
    APqua = pqua(t[idx1:idx2])
    APtout = t[idx0:idx2]
    #print(APlin)
    #APout = np.concatenate(APlin[:-1],APqua[:])
    APout = np.array(list(APlin[:])+list(APqua[:]))
    return APout, APtout
    


# In[5]:

pathToFiles = "./GroupI/"  # change path to files here
simFileNames = "Well*Recording*.csv"
SimFiles = glob.glob(pathToFiles+simFileNames)
SimFiles.remove("./GroupI/Well2Recording3.csv")
percentile = 0.95
x = 0.8


# In[48]:

get_ipython().magic('matplotlib nbagg')


# In[58]:

for simFile in SimFiles:
    AP = np.loadtxt(simFile, delimiter=',')
    numIdx = int(percentile*len(AP))
    maxSig = np.sort(AP[:,1])[numIdx]
    AP[:,1] = AP[:,1]/maxSig
    (a,b)=getSmoothAP(AP[:,1],AP[:,0],t0=0.047,t1=0.054,t2=0.35)
    #APDx = getAPDx(AP[:,1],AP[:,0],check2groups=False)
    APDx = getAPDx(a,b,check2groups=False)
    if APDx[0] != np.inf:
        APar = getAPar(a,b,APDx[0],APDx[2],APDx[4])
    else:
        APar = np.inf
    print(APDx[0],APar)
    #plt.plot(AP[:,0],AP[:,1])
    #plt.plot(b,a)
    plt.plot(APDx[0]*1000.,APar,'ro')
plt.title('Group I optical mapping AP data')
plt.xlabel('APD [ms]')
plt.ylabel('AP Area Ratio')
plt.show()


# In[57]:

# test getsmoothAP()
AP = np.loadtxt(SimFiles[0], delimiter=',')
print( np.where(AP[:,0]==0.075)[0])
(a,b)=getSmoothAP(AP[:,1],AP[:,0],t0=0.047,t1=0.054,t2=0.35)
plt.plot(AP[:,0]*1000.,AP[:,1],label='raw')
plt.plot(b*1000.,a,label='filtered')
plt.title('Group I optical mapping AP data filtered by fitting')
plt.xlabel('t [ms]')
plt.ylabel('Vm [AU]')
plt.legend()
plt.show()

