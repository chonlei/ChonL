import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

simFilesDir = "./outputPaciSimulation*.txt"  # change path to files here
SimFiles = glob.glob(simFilesDir)
orginalModel = "outputPaciSimulation1-32.txt"

print "Plotting simulation results..."
for simFile in SimFiles:
    #print simFile
    Vt = np.loadtxt(simFile)
    plt.plot(Vt[:,0]*10., Vt[:,1],alpha=0.5)#,label=simFile)
    if orginalModel in simFile:
        plt.plot(Vt[:,0]*10., Vt[:,1],'k',linewidth=2,label="Original Model")
        print "Found original model and plotted"


print "Finished plotting."
plt.xlabel("t [ms]")
plt.ylabel("V [mV]")
plt.legend()
if len(sys.argv)==2 and sys.argv[1]=="--save":
    savefigname = "outputPaciSimulation.png"
    plt.savefig(savefigname)
    print "Saved figure as %s"%savefigname
plt.show()

