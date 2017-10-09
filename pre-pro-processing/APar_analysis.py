import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

pathToFiles = "./out_PaciSweep_1_corrected/"  # change path to files here
simFileNames = "outputPaciSimulation*.txt"
SimFiles = glob.glob(pathToFiles+simFileNames)
percentile = 0.95
x = 0.8

sweepFileName = "conductanceDataTest*.txt"
sweeping_param = np.loadtxt(glob.glob(pathToFiles+sweepFileName)[0])
paramIdx = {1:r"$g_{Na}$",2:r"$g_{CaL}$",3:r"$g_{Kr}$"}
paramIdxTxt = {1:"gNa",2:"gCaL",3:"gKr"}


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


###############################################################################
'''
APDX = []
for simFile in SimFiles:
    #print simFile
    AP = np.loadtxt(simFile)
    APDx = getAPDx(AP[:,1],AP[:,0],percentile,x)
    APDX.append(APDx[0])

APDX = np.array(APDX)
plt.hist(APDX[APDX!=np.inf])
'''

###############################################################################
## Plot 1
# AP for varying only one parameter at a time

## Na
plt.figure(1)
listOfIdx = []
for i in range(len(sweeping_param)):
    if sweeping_param[i,2] == sweeping_param[i,3] == 1:
        listOfIdx.append(i)

for i in listOfIdx:
    apfile = [f for f in glob.glob(pathToFiles+"outputPaciSimulation%d*"%i) if "Simulation%d.txt"%i in f or "Simulation%d_nonSS.txt"%i in f][0]
    AP = np.loadtxt(apfile)
    plt.plot(AP[:,0],AP[:,1],label=r'%%$g_{Na}=%g$'%sweeping_param[i,1])
plt.legend()
plt.xlabel("t [ms]")
plt.ylabel("V [mV]")
plt.title(r"Vary solely $g_{Na}$ by indicated ratio")
plt.ylim([-85, 45])
plt.savefig("vgNa_only.png")

## CaL
plt.figure(2)
listOfIdx = []
for i in range(len(sweeping_param)):
    if sweeping_param[i,1] == sweeping_param[i,3] == 1:
        listOfIdx.append(i)

for i in listOfIdx:
    apfile = [f for f in glob.glob(pathToFiles+"outputPaciSimulation%d*"%i) if "Simulation%d.txt"%i in f or "Simulation%d_nonSS.txt"%i in f][0]
    AP = np.loadtxt(apfile)
    plt.plot(AP[:,0],AP[:,1],label=r"%%$g_{CaL}=%g$"%sweeping_param[i,2])
plt.legend()
plt.xlabel("t [ms]")
plt.ylabel("V [mV]")
plt.title(r"Vary solely $g_{CaL}$ by indicated ratio")
plt.ylim([-85, 45])
plt.savefig("vgCaL_only.png")

## Kr
plt.figure(3)
listOfIdx = []
for i in range(len(sweeping_param)):
    if sweeping_param[i,1] == sweeping_param[i,2] == 1:
        listOfIdx.append(i)

for i in listOfIdx:
    apfile = [f for f in glob.glob(pathToFiles+"outputPaciSimulation%d*"%i) if "Simulation%d.txt"%i in f or "Simulation%d_nonSS.txt"%i in f][0]
    AP = np.loadtxt(apfile)
    plt.plot(AP[:,0],AP[:,1],label=r"%%$g_{Kr}=%g$"%sweeping_param[i,3])
plt.legend()
plt.xlabel("t [ms]")
plt.ylabel("V [mV]")
plt.title(r"Vary solely $g_{Kr}$ by indicated ratio")
plt.ylim([-85, 45])
plt.savefig("vgKr_only.png")


###############################################################################
## Plot 2
# APD90 for varying only one parameter at a time

plt.figure(4)
scanRange = [0.1,0.5,1.0,1.5,2.0,2.5,3.0]
pltStyle = ['bo','g^','rs','cv','mX','kP','yD']
# params idx, e.g.[1,2,3] for 3 params indices (0 reserved for combination idx)
paramIdxx = np.arange(1,len(paramIdx)+1)

for rr in xrange(len(scanRange)):
    for i in xrange(len(sweeping_param)):
        for j in paramIdxx:
            othersAre1 = True
            for k in paramIdxx[paramIdxx!=j]:
                othersAre1 = (othersAre1 and sweeping_param[i,k]==1.0)
            if sweeping_param[i,j] == scanRange[rr] and othersAre1:
                apfile = [f for f in glob.glob(pathToFiles+"outputPaciSimulation%d*"%i) if "Simulation%d.txt"%i in f or "Simulation%d_nonSS.txt"%i in f][0]
                AP = np.loadtxt(apfile)
                APDx = getAPDx(AP[:,1],AP[:,0])
                if APDx[0] != np.inf:
                    APar = getAPar(AP[:,1],AP[:,0],APDx[0],APDx[2],APDx[4])
                else:
                    APar = np.inf
                
                # have label once only
                if j == 1:
                    plt.plot(j,APar,pltStyle[rr],ms=10,label="%%=%g"%scanRange[rr])
                else:
                    plt.plot(j,APar,pltStyle[rr],ms=10)

plt.legend()
plt.ylabel(r"AP area ratio")
plt.xlim([0.8,4])
plt.xticks(paramIdx.keys(), paramIdx.values())
plt.title(r"AP area ratio defined by $APD_{90}$ for varying only one parameter at a time")
plt.savefig("APD90-def-APar_1at1time.png")


###############################################################################
## Plot 3
# APD90 for varying only two parameters at a time
plt.figure(5)

## change gNa , gCaL
holdConst = (3,1.0)  # (paramIdx,valueToHoldAt)
toVary = (1,2)

## change gNa , gCaL
holdConst = (2,1.0)  # (paramIdx,valueToHoldAt)
toVary = (1,3)

## change gNa , gCaL
holdConst = (1,1.0)  # (paramIdx,valueToHoldAt)
toVary = (2,3)

Idxmat = np.zeros((len(scanRange),len(scanRange)))
Aramat = np.zeros((len(scanRange),len(scanRange)))

for i in xrange(len(sweeping_param)):
    if sweeping_param[i,holdConst[0]] == holdConst[1]:
        ii = scanRange.index(sweeping_param[i,toVary[0]])
        jj = scanRange.index(sweeping_param[i,toVary[1]])
        Idxmat[ii,jj] = sweeping_param[i,0]  #just get idx
        
        apfile = [f for f in glob.glob(pathToFiles+"outputPaciSimulation%d*"%i) if "Simulation%d.txt"%i in f or "Simulation%d_nonSS.txt"%i in f][0]
        #print i, apfile
        AP = np.loadtxt(apfile)
        APDx = getAPDx(AP[:,1],AP[:,0])
        if APDx[0] != np.inf and i!=18: #idx 18 here is a bit weird
            print i, APDx[0],APDx[2],APDx[4], AP[0,1]
            Aramat[ii,jj] = getAPar(AP[:,1],AP[:,0],APDx[0],APDx[2],APDx[4])
        else:
            Aramat[ii,jj] = np.inf

plt.imshow(Aramat.transpose(),cmap='hot', interpolation='nearest')  #,vmax=0.8,vmin=0.1
plt.xticks(range(len(scanRange)), scanRange)
plt.yticks(range(len(scanRange)), scanRange)
plt.xlabel("%"+paramIdx[toVary[0]])
plt.ylabel("%"+paramIdx[toVary[1]])
plt.title(r"AP area ratio defined by $APD_{90}$ [ms] as a function of two parameters")
plt.colorbar()
plt.savefig("APD90-def-APar_%s_%s.png"%(paramIdxTxt[toVary[0]],paramIdxTxt[toVary[1]]))

# just to verify manually the V-t trace
i = Idxmat[0,1]
apfile = [f for f in glob.glob(pathToFiles+"outputPaciSimulation%d*"%i) if "Simulation%d.txt"%i in f or "Simulation%d_nonSS.txt"%i in f][0]
AP = np.loadtxt(apfile)
plt.plot(AP[:,0],AP[:,1])


