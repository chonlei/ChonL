import re
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

for loop_idx in [0]:#[0,1,2,3,11,12,13,21,22,23,31,32]:
	pathToSave = "./-100-600analysis/"
	saveName = "IKr_%d"%loop_idx
	#niceName = "I_{K}"
	#It_V_idx = 7
	#loopList = [("IKr","I_{Kr}",1),("IKs","I_{Ks}",2),("IKto","I_{Kto}",3),("","",4),("","",5),("","",6)]
        loopList = [("IKs","I_{Ks}",2)]
	isGetMin = False
	checkIsClamped = [False, 1,2,3]
	pathToFiles = "./tempout/"
	pCSimFiles = "loop%doutPaciICaL*.txt"%loop_idx
	pCSimReV = ".*loop%doutPaciICaL(.+)\.txt"%loop_idx
	timeWindow = (-np.inf,np.inf)

	PCSims = glob.glob(pathToFiles+pCSimFiles)
	PCSims.sort()
	t_idx = 0
	dt = 1

	V = []
	I = []

	plt.figure(2)
	#plt.title(r"$%s$ excluding $I_{Na}$ and $I_{CaL}$ patch clamp raw traces"%"I_{all}")
	plt.xlabel("$t$ [ms]")
	plt.ylabel(r"$I$ [A/F]")

	for pCSim in PCSims:
	    # extract V from file name
	    V.append(float(re.findall(pCSimReV,pCSim)[0]))
	    # load simulation
	    It_V = np.loadtxt(pCSim)
	    # define window of interest
	    tw = ((It_V[:,0]<timeWindow[1]) & (It_V[:,0]>timeWindow[0]))
            #print tw
	    # get indices of currents of interest
	    IdxIint = [loopList[i][2] for i in range(len(loopList))]
	    # calculate the sum of all interested currents
	    Isum = np.sum(It_V[:,IdxIint],1)
            #print Isum
	    # user defined to find min or max of the current (depends on type of current, i.e. inward or outward current)
	    if isGetMin:
		indexOfInt = np.argmin(Isum[tw])
	    else:
		indexOfInt = np.argmax(Isum[tw])
	    I.append(It_V[tw,:][:,IdxIint][indexOfInt,:])
	    plt.plot(It_V[:,t_idx][tw]*(1./dt),Isum[tw],label="V held@%g mV"%V[-1])
	    # check clamped values vary less than 1e-10
	    if checkIsClamped[0]:
		for checkClampIdx in checkIsClamped[1:]:
		    if np.abs(np.max(It_V[:,checkClampIdx])-np.min(It_V[:,checkClampIdx])) > 1e-10:
			print "WARNING: file %s has clamp value (column index=%d) varies more than 1e-10"%(pCSim,checkClampIdx)
        #plt.ylim([-1,15])
	plt.savefig("%s%s_PC_raw_sim.png"%(pathToSave,saveName))

        '''
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
	np.random.seed(97)
	def get_spaced_colors(n):
	    max_value = 16581375 #255**3
	    interval = int(max_value / n)
	    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]
	    max_value = float(255.)
	    return [(int(i[:2], 16)/max_value, int(i[2:4], 16)/max_value, int(i[4:], 16)/max_value) for i in colors]
	stack_coll = plt.stackplot(Vsorted,np.array(Isorted).transpose(), colors=[np.random.rand(3,1) for i in range(len(loopList))])
	#stack_coll = plt.stackplot(Vsorted,np.array(Isorted).transpose(), colors=get_spaced_colors(len(loopList)))

	# make legend
	proxy_rects = [Rectangle((0, 0), 1, 1, fc=pc.get_facecolor()[0]) for pc in stack_coll]
	plt.legend(proxy_rects, [r"$%s$"%loopList[i][1] for i in range(len(loopList))], loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3, fancybox=True, shadow=True)

	#plt.show()
	#plt.legend()
	plt.savefig("%s%s_IV_curve.png"%(pathToSave,saveName))

	np.savetxt("%s%s_I.txt"%(pathToSave,saveName),np.array(Isorted).transpose())
	np.savetxt("%s%s_V.txt"%(pathToSave,saveName),np.array(Vsorted))'''
	plt.close('all')


