import re
import glob
import numpy as np
import matplotlib.pyplot as plt

pathToSave = "./outdata/ICaLIVexp/"
saveName = "IK"
niceName = "I_{K}"
It_V_idx = 7
#loopList = [("IK1","I_{K1}",1),("IKr","I_{Kr}",2),("IKs","I_{Ks}",3),("IKto","I_{Kto}",4),("ICap","I_{Cap}",5),("If","I_{f}",6),("INaCa","I_{NaCa}",7),("INaK","I_{NaK}",8),("ICaL","I_{CaL}",9),("INa","I_{Na}",10)]
#loopList = [("INa","I_{Na}",4)]
loopList = [("ICaL","I_{CaL}",2)]
isGetMin = True
checkIsClamped = [False, 1,2,3]
pathToFiles = "./outdata/ICaLIVexp/" #"/home/scratch/chaste/testoutput/out_PaciSweep/" #"./outdata/INaIVexp/" #"./outdata/ICaLIVexp/"
pCSimFiles = "outPaciICaL*.txt"
pCSimReV = ".*outPaciICaL(.+)\.txt"
timeWindow = (-np.inf,np.inf)

PCSims = glob.glob(pathToFiles+pCSimFiles)
PCSims.sort()
t_idx = 0
dt = 1

for (saveName,niceName,It_V_idx) in loopList:
	V = []
	I = []

	plt.figure(1)
	plt.title(r"$%s$ patch clamp raw traces"%niceName)
	plt.xlabel("t [ms]")
	plt.ylabel(r"$%s$ [A/F]"%niceName)
	for pCSim in PCSims:
	    # extract V from file name
	    V.append(float(re.findall(pCSimReV,pCSim)[0]))
	    # load simulation
	    It_V = np.loadtxt(pCSim)
	    # user defined to find min or max of the current (depends on type of current, i.e. inward or outward current)
	    if isGetMin:
		I.append(np.min(It_V[:,It_V_idx]))
	    else:
		I.append(np.max(It_V[:,It_V_idx][((It_V[:,0]<timeWindow[1]) & (It_V[:,0]>timeWindow[0]) )]))
	    # check clamped values vary less than 1e-10
	    if checkIsClamped[0]:
		for checkClampIdx in checkIsClamped[1:]:
		    if np.abs(np.max(It_V[:,checkClampIdx])-np.min(It_V[:,checkClampIdx])) > 1e-10:
		        print "WARNING: file %s has clamp value (column index=%d) varies more than 1e-10"%(pCSim,checkClampIdx)
	    plt.plot(It_V[:,t_idx][((It_V[:,0]<timeWindow[1]) & (It_V[:,0]>timeWindow[0]) )]*(1./dt),It_V[:,It_V_idx][((It_V[:,0]<timeWindow[1]) & (It_V[:,0]>timeWindow[0]) )],label="V held@%g mV"%V[-1])
	#plt.legend()
	plt.savefig("%s%s_PC_raw_sim.png"%(pathToSave,saveName))

	# sort in the order of V
	Isorted = [i for (v,i) in sorted(zip(V,I))]
	Vsorted = sorted(V)

	plt.figure(2)
	plt.title(r"$%s$ IV curve"%niceName)
	plt.xlabel("V [mV]")
	plt.ylabel(r"$sign(%s)max(|%s}|)$ [A/F]"%(niceName,niceName))
	plt.plot(Vsorted,Isorted,'bo-')
	plt.savefig("%s%s_IV_curve.png"%(pathToSave,saveName))
	np.savetxt("%s%s_IV_curve.txt"%(pathToSave,saveName),np.array([V,I]))
        plt.close('all')


