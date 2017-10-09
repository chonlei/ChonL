import re
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

pathToSave = "./allCurrents/"
saveName = "stack_all_K_wholeStep"
#niceName = "I_{K}"
#It_V_idx = 7
#loopList = [("IK1","I_{K1}",1,False),("IKr","I_{Kr}",2,False),("IKs","I_{Ks}",3,False),("IKto","I_{Kto}",4,False),("ICap","I_{Cap}",5,False),("If","I_{f}",6,False),("INaCa","I_{NaCa}",7,False),("INaK","I_{NaK}",8,False),("ICaL","I_{CaL}",9,True),("INa","I_{Na}",10,True)]
loopList = [("IK1","I_{K1}",1,False),("IKr","I_{Kr}",2,False),("IKs","I_{Ks}",3,False),("IKto","I_{Kto}",4,False)]
#loopList = [("IK1","I_{K1}",1,False),("IKr","I_{Kr}",2,False),("IKs","I_{Ks}",3,False),("IKto","I_{Kto}",4,False),("ICap","I_{Cap}",5,False),("If","I_{f}",6,False),("INaCa","I_{NaCa}",7,False),("INaK","I_{NaK}",8,False)]
#loopList = [("IK1","I_{K1}",1),("IKr","I_{Kr}",2),("IKs","I_{Ks}",3),("IKto","I_{Kto}",4)]
#loopList = [("IK1","I_{K1}",1),("IKr","I_{Kr}",2),("IKs","I_{Ks}",3),("IKto","I_{Kto}",4),("ICap","I_{Cap}",5),("If","I_{f}",6),("INaCa","I_{NaCa}",7),("INaK","I_{NaK}",8)]
isGetMin = False
checkIsClamped = [False, 1,2,3]
pathToFiles = "./tempout/"
pCSimFiles = "outPaciICaL*.txt"
pCSimReV = ".*outPaciICaL(.+)\.txt"
timeWindow = (0,500)

PCSims = glob.glob(pathToFiles+pCSimFiles)
PCSims.sort()
t_idx = 0
dt = 1

V = []
I = []

for pCSim in PCSims:
    # extract V from file name
    V.append(float(re.findall(pCSimReV,pCSim)[0]))
    # load simulation
    It_V = np.loadtxt(pCSim)
    # define window of interest
    tw = ((It_V[:,0]<timeWindow[1]) & (It_V[:,0]>timeWindow[0]))
    # get indices of currents of interest
    IdxIint = [loopList[i][2] for i in range(len(loopList))]
    # calculate the sum of all interested currents
    #Isum = np.sum(It_V[:,IdxIint],1)
    # user defined to find min or max of the current (depends on type of current, i.e. inward or outward current)
    I.append([])
    for i in xrange(len(IdxIint)):
        if loopList[i][3]:
            indexOfInt = np.argmin(It_V[tw,:][:,IdxIint[i]])
        else:
            indexOfInt = np.argmax(It_V[tw,:][:,IdxIint[i]])
        I[-1].append(It_V[tw,:][indexOfInt,IdxIint[i]])
    # check clamped values vary less than 1e-10
    if checkIsClamped[0]:
	for checkClampIdx in checkIsClamped[1:]:
	    if np.abs(np.max(It_V[:,checkClampIdx])-np.min(It_V[:,checkClampIdx])) > 1e-10:
	        print "WARNING: file %s has clamp value (column index=%d) varies more than 1e-10"%(pCSim,checkClampIdx)

# sort in the order of V
Isorted = [i for (v,i) in sorted(zip(V,I))]
Vsorted = sorted(V)

plt.figure(1)
#plt.title(r"Stackplot for currents contribution in IV curve")
plt.xlabel("V [mV]")
plt.ylabel(r"I [A/F]")
# eqv. to area plot in Matlab
#stack_coll = plt.stackplot(Vsorted,np.array(Isorted).transpose())
# set colours
#np.random.seed(97)
def get_spaced_colors(n):
    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]
    max_value = float(255.)
    return [(int(i[:2], 16)/max_value, int(i[2:4], 16)/max_value, int(i[4:], 16)/max_value) for i in colors]
#stack_coll = plt.stackplot(Vsorted,np.array(Isorted).transpose(), colors=[np.random.rand(3,1) for i in range(len(loopList))])
stack_coll = plt.stackplot(Vsorted,np.array(Isorted).transpose(), colors=get_spaced_colors(len(loopList)))

# make legend
proxy_rects = [Rectangle((0, 0), 1, 1, fc=pc.get_facecolor()[0]) for pc in stack_coll]
plt.legend(proxy_rects, [r"$%s$"%loopList[i][1] for i in range(len(loopList))], loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3, fancybox=True, shadow=True)

#plt.show()
#plt.legend()
plt.savefig("%s%s_IV_curve.png"%(pathToSave,saveName))

np.savetxt("%s%s_I.txt"%(pathToSave,saveName),np.array(Isorted).transpose())
np.savetxt("%s%s_V.txt"%(pathToSave,saveName),np.array(Vsorted))
plt.close('all')


