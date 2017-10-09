import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import re
from scipy.optimize import basinhopping


#pathToData = "./leak_corrected_data/"
#dataFile = "corrected_Tr_140202_Chip_0_E3_NaIV.txt"
pathToData = "./leak_corrected_data_IK/"
dataFile = "corrected_Tr_140202_Chip_0_B3_KIV_IKsr.txt"

pathToFiles = "./sim_Kr_Ks_tau03_Nai10/"
pCSimFiles = "*outPaciICaL*.txt"
pCSimReV = ".*outPaciICaL(.+)\.txt"
timeWindow = (400,500)  # ms


PCSims = glob.glob(pathToFiles+pCSimFiles)
PCSims.sort()
#loopList = [("IKr","I_{Kr}",1),("IKs","I_{Ks}",2),("IKto","I_{Kto}",3),("NaCa","",4),("NaK","",5),("INa","",6)]
t_idx = 0
kr_idx = 1
ks_idx = 2
na_idx = 6
dt = 1

voltageFile = "./sweep_voltage.txt" # change path to files here
dataST = 50  # [ms] data step potential starting time




def get_spaced_colors(n):
    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]
    max_value = float(255.)
    return [(int(i[:2], 16)/max_value, int(i[2:4], 16)/max_value, int(i[4:], 16)/max_value) for i in colors]





def getMaxCurrent0(PCSims, dataFile, pathToData, fixVat50=False):
    # Return maximum current for IKr and IKs in data and simulation
    data = np.loadtxt(pathToData+dataFile)

    dataV = np.loadtxt(voltageFile)
    # assumed data time unit [us]
    dataTW = ( (data[:,0]>(timeWindow[0]+dataST)*1e3) & (data[:,0]<(timeWindow[1]+dataST)*1e3) )
    
    maxIdata = 0
    maxIKrsim = 0
    maxIKssim = 0
    
    for pCSim in PCSims:
        Vsim = float(re.findall(pCSimReV,pCSim)[0]) - 0.01
        #if Vsim not in dataV:
        if fixVat50:
            if Vsim != 50.0:
                continue
        else:
            if (Vsim not in dataV):
                continue
        
        # else do something
        It = np.loadtxt(pCSim)
        data_idx = list(dataV).index(Vsim) + 1  # 0 is time
        
        ItTW = ( (It[:,0]>timeWindow[0]) & (It[:,0]<timeWindow[1]) )
        
        tempId = np.max(data[dataTW,data_idx])
        tempIKrs = np.max(It[ItTW,kr_idx])
        tempIKss = np.max(It[ItTW,ks_idx])
        
        maxIdata = tempId if (tempId>maxIdata) else maxIdata
        maxIKrsim = tempIKrs if (tempIKrs>maxIKrsim) else maxIKrsim
        maxIKssim = tempIKss if (tempIKss>maxIKssim) else maxIKssim
    return maxIdata,maxIKrsim,maxIKssim





def getMaxCurrent(PCSims, dataFile, pathToData, scanningRrKrKs=1.0, absR=None, fixVat50=False):
    # Return maximum current for IKr and IKs in data and simulation given the ratio of IKr and IKs
    if absR == None:
        absRtemp = getMaxCurrent0(PCSims, dataFile, pathToData, fixVat50)
        absR = (absRtemp[0]/absRtemp[1], absRtemp[0]/absRtemp[2])
        #print "Using (RKr,RKs) = ", absR
    if scanningRrKrKs>1.0 or scanningRrKrKs<0:
        raise Exception("ratio not within range [0,1].")
    
    
    data = np.loadtxt(pathToData+dataFile)

    dataV = np.loadtxt(voltageFile)
    # assumed data time unit [us]
    dataTW = ( (data[:,0]>(timeWindow[0]+dataST)*1e3) & (data[:,0]<(timeWindow[1]+dataST)*1e3) )
    
    maxIdata = 0
    maxIsim = 0
    
    for pCSim in PCSims:
        Vsim = float(re.findall(pCSimReV,pCSim)[0]) - 0.01
        #if Vsim not in dataV:
        if fixVat50:
            if Vsim != 50.0:
                continue
        else:
            if (Vsim not in dataV):
                continue
        
        # else do something
        It = np.loadtxt(pCSim)
        data_idx = list(dataV).index(Vsim) + 1  # 0 is time
        
        ItTW = ( (It[:,0]>timeWindow[0]) & (It[:,0]<timeWindow[1]) )
        
        tempId = np.max(data[dataTW,data_idx])
        
        IKrSIM = It[ItTW,kr_idx]*absR[0]*scanningRrKrKs
        IKsSIM = It[ItTW,ks_idx]*absR[1]*(1.0-scanningRrKrKs)
        tempIs = np.max(IKrSIM + IKsSIM)
        
        maxIdata = tempId if (tempId>maxIdata) else maxIdata
        maxIsim = tempIs if (tempIs>maxIsim) else maxIsim
    
    return maxIdata,maxIsim





def getCurrentTraces(PCSims, dataFile, pathToData, scanningRrKrKs=1.0, fixVat50=False, save=None):
    absRtemp = getMaxCurrent0(PCSims, dataFile, pathToData, fixVat50)
    absR = (absRtemp[0]/absRtemp[1], absRtemp[0]/absRtemp[2])
    maxC = getMaxCurrent(PCSims, dataFile, pathToData, rr, absR=absR, fixVat50=fixVat50)
    
    dataV = np.loadtxt(voltageFile)
    colors = get_spaced_colors(len(dataV))
    ccount = 0
    for pCSim in PCSims:
        Vsim = float(re.findall(pCSimReV,pCSim)[0]) - 0.01
        if (Vsim not in dataV):
            continue
        # else do something
        It = np.loadtxt(pCSim)
        IKrSIM = It[:,kr_idx]*absR[0]*scanningRrKrKs * maxC[0]/maxC[1]
        IKsSIM = It[:,ks_idx]*absR[1]*(1.0-scanningRrKrKs) * maxC[0]/maxC[1]
        plt.plot(It[:,0],IKrSIM+IKsSIM,color=colors[ccount],label='V=%gmV'%Vsim)
        ccount += 1
    plt.xlabel("t [ms]")
    plt.ylabel("I [A/F]")
    plt.ylim([-1, absRtemp[0]*1.4])
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3, fancybox=True, shadow=True)
    if save==None:
        plt.show()
    else:
        plt.savefig("perIKr=%g_fixV50=%s_%s.png"%(scanningRrKrKs,fixVat50,save))
    plt.close('all')





def compareNaTraces(PCSims, dataFile, pathToData,save=None):
    data = np.loadtxt(pathToData+dataFile)
    dataV = np.loadtxt(voltageFile)
    colors = get_spaced_colors(len(dataV))
    ccount = 0
    for pCSim in PCSims:
        Vsim = float(re.findall(pCSimReV,pCSim)[0]) - 0.01
        if (Vsim not in dataV) or Vsim==-30.0 or Vsim==-10.0 or Vsim==10.0 or Vsim==30.0:
            continue
        # else do something
        It = np.loadtxt(pCSim)
        data_idx = list(dataV).index(Vsim) + 1  # 0 is time
        # It[:,na_idx]+It[:,3]+It[:,1]*20+It[:,2]*10+It[:,4]+It[:,5];   It[:,3]*18;   It[:,1]*30
        ### It[:,3]*8+It[:,na_idx]+It[:,2]*10+It[:,1]*10   +   It[:,4]*4+It[:,5]*5
        plt.plot(It[:,0],It[:,3]*7+It[:,na_idx]+It[:,2]*10   +   It[:,4]*2,color=colors[ccount],linewidth=2.0)#,label='V=%gmV'%Vsim)
        plt.plot(data[:,0]*1e-3-50, data[:,data_idx],color=colors[ccount])#,label='V=%gmV'%Vsim)
        ccount += 1
    plt.xlabel("t [ms]")
    plt.ylabel("I [A/F]")
    plt.xlim([-30, 60])
    plt.ylim([-10, 15])
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3, fancybox=True, shadow=True)
    if save==None:
        plt.show()
    else:
        plt.savefig("perIKr=%g_fixV50=%s_%s.png"%(scanningRrKrKs,fixVat50,save))
    plt.close('all')





def compareNaTracesDiff(PCSims, dataFile, pathToData,save=None):
    data = np.loadtxt(pathToData+dataFile)
    dataV = np.loadtxt(voltageFile)
    colors = get_spaced_colors(len(dataV))
    ccount = 0
    for pCSim in PCSims:
        Vsim = float(re.findall(pCSimReV,pCSim)[0]) - 0.01
        if (Vsim not in dataV) or Vsim==-30.0 or Vsim==-10.0 or Vsim==10.0 or Vsim==30.0:# or Vsim!=40.0:
            continue
        # else do something
        It = np.loadtxt(pCSim)
        data_idx = list(dataV).index(Vsim) + 1  # 0 is time
        
        #### match the time
        ItTW = ( (It[:,0]>0.5) & (It[:,0]<10) )
        timeData = [round(i*1e-3-50,3) for i in data[:,0]]
        ItTWData = [i for i, item in enumerate(It[ItTW,0]) if round(item,3) in set(timeData)]
        It = It[ItTW,:][ItTWData,:]
        
        dataTWIt = [i for i, item in enumerate(data[:,0]*1e-3-50) if round(item,3) in set(It[:,0])]
        data = data[dataTWIt,:]
        
        # define the objective function
        def f_iNa_g_tomin(x):
            (rto,rkr,rncx,rnak,rna) = x
            I1 = It[:,1]*rkr + It[:,3]*rto + It[:,4]*rncx + It[:,5]*rnak + It[:,6]*rna
            I2 = data[:,data_idx]
            return np.linalg.norm(I1-I2)
        
        # the starting point
        x0 = [1., 7., 2., 1., 1.]
        
        # the bounds
        xmin = [0., 0., 0., 0., 0.] #[ 0., 18.55682001,  0. ,  0. , 0.78500075]#
        xmax = [100., 20., 11., 50., 11.]
        
        # rewrite the bounds in the way required by L-BFGS-B
        bounds = [(low, high) for low, high in zip(xmin, xmax)]
        
        # use method L-BFGS-B because the problem is smooth and bounded
        minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds)
        res = basinhopping(f_iNa_g_tomin, x0, minimizer_kwargs=minimizer_kwargs)
        
        print "global minimum: x = ", res.x, ", f(x0) = ", res.fun
        (rtom,rkrm,rncxm,rnakm,rnam) = res.x
        I1m = It[:,1]*rkrm + It[:,3]*rtom + It[:,4]*rncxm + It[:,5]*rnakm + It[:,6]*rnam
        
        
        # It[:,na_idx]+It[:,3]+It[:,1]*20+It[:,2]*10+It[:,4]+It[:,5];   It[:,3]*18;   It[:,1]*30
        ### It[:,3]*8+It[:,na_idx]+It[:,2]*10+It[:,1]*10   +   It[:,4]*4+It[:,5]*5
        #plt.plot(It[:,0],It[:,3]*7+It[:,na_idx]+It[:,2]*10   +   It[:,4]*2,color=colors[ccount],linewidth=2.0)#,label='V=%gmV'%Vsim)
        plt.plot(It[:,0], I1m, color=colors[ccount])
        plt.plot(data[:,0]*1e-3-50, data[:,data_idx],color=colors[ccount])#,label='V=%gmV'%Vsim)
        ccount += 1
    plt.xlabel("t [ms]")
    plt.ylabel("I [A/F]")
    plt.xlim([-30, 60])
    plt.ylim([-10, 15])
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3, fancybox=True, shadow=True)
    if save==None:
        plt.show()
    else:
        plt.savefig("perIKr=%g_fixV50=%s_%s.png"%(scanningRrKrKs,fixVat50,save))
    plt.close('all')






def compareNaTracesDiff2(PCSims, dataFile, pathToData,save=None):
    ## Compare all voltage
    # load data for minimisation and format correction
    data = np.loadtxt(pathToData+dataFile)
    dataV = np.loadtxt(voltageFile)
    AllIt = []
    AllData = []
    for pCSim in PCSims:
        Vsim = float(re.findall(pCSimReV,pCSim)[0]) - 0.01
        if (Vsim not in dataV): # or Vsim==-30.0 or Vsim==-10.0 or Vsim==10.0 or Vsim==30.0:# or Vsim!=40.0:
            continue
        # else do something
        It = np.loadtxt(pCSim)
        data_idx = list(dataV).index(Vsim) + 1  # 0 is time
        #### match the time
        ItTW = ( (It[:,0]>0.5) & (It[:,0]<10) )
        timeData = [round(i*1e-3-50,3) for i in data[:,0]]
        ItTWData = [i for i, item in enumerate(It[ItTW,0]) if round(item,3) in set(timeData)]
        It = It[ItTW,:][ItTWData,:]
        
        dataTWIt = [i for i, item in enumerate(data[:,0]*1e-3-50) if round(item,3) in set(It[:,0])]
        data = data[dataTWIt,:]
        AllIt.append(It)
        AllData.append(data[:,data_idx])

    # define the objective function
    def f_iNa_g_tomin(x):
        (rto,rkr,rncx,rnak,rna) = x
        ERR = 0.0
        for i in range(len(AllData)):
            I1 = AllIt[i][:,1]*rkr + AllIt[i][:,3]*rto + AllIt[i][:,4]*rncx + AllIt[i][:,5]*rnak + AllIt[i][:,6]*rna
            I2 = AllData[i][:]
            ERR += np.linalg.norm(I1-I2)
        return ERR
    
    # the starting point
    x0 = [1., 7., 2., 1., 1.]
    
    # the bounds
    xmin = [0., 0., 0., 0., 0.] #[ 0., 18.55682001,  0. ,  0. , 0.78500075]#
    xmax = [100., 20., 11., 50., 11.]
    
    # rewrite the bounds in the way required by L-BFGS-B
    bounds = [(low, high) for low, high in zip(xmin, xmax)]
    
    # use method L-BFGS-B because the problem is smooth and bounded
    minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds)
    res = basinhopping(f_iNa_g_tomin, x0, minimizer_kwargs=minimizer_kwargs)
    
    print "global minimum: x = ", res.x, ", f(x0) = ", res.fun
    (rtom,rkrm,rncxm,rnakm,rnam) = res.x
    # Reload data for plotting
    data = np.loadtxt(pathToData+dataFile)
    dataV = np.loadtxt(voltageFile)
    colors = get_spaced_colors(len(dataV))
    ccount = 0
    for pCSim in PCSims:
        Vsim = float(re.findall(pCSimReV,pCSim)[0]) - 0.01
        if (Vsim not in dataV) or Vsim==-30.0 or Vsim==-10.0 or Vsim==10.0 or Vsim==30.0:# or Vsim!=40.0:
            continue
        # else do something
        It = np.loadtxt(pCSim)
        data_idx = list(dataV).index(Vsim) + 1  # 0 is time
        I1m = It[:,1]*rkrm + It[:,3]*rtom + It[:,4]*rncxm + It[:,5]*rnakm + It[:,6]*rnam
        
        
        # It[:,na_idx]+It[:,3]+It[:,1]*20+It[:,2]*10+It[:,4]+It[:,5];   It[:,3]*18;   It[:,1]*30
        ### It[:,3]*8+It[:,na_idx]+It[:,2]*10+It[:,1]*10   +   It[:,4]*4+It[:,5]*5
        #plt.plot(It[:,0],It[:,3]*7+It[:,na_idx]+It[:,2]*10   +   It[:,4]*2,color=colors[ccount],linewidth=2.0)#,label='V=%gmV'%Vsim)
        plt.plot(It[:,0], I1m, color=colors[ccount])
        plt.plot(data[:,0]*1e-3-50, data[:,data_idx],color=colors[ccount])#,label='V=%gmV'%Vsim)
        ccount += 1
    plt.xlabel("t [ms]")
    plt.ylabel("I [A/F]")
    plt.xlim([-30, 60])
    plt.ylim([-10, 15])
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3, fancybox=True, shadow=True)
    if save==None:
        plt.show()
    else:
        plt.savefig("perIKr=%g_fixV50=%s_%s.png"%(scanningRrKrKs,fixVat50,save))
    plt.close('all')








def compareKTracesDiff2(PCSims, dataFile, pathToData,save=None):
    ## Compare all voltage
    # load data for minimisation and format correction
    data = np.loadtxt(pathToData+dataFile)
    dataV = np.loadtxt(voltageFile)
    AllIt = []
    AllData = []
    for pCSim in PCSims:
        Vsim = float(re.findall(pCSimReV,pCSim)[0]) - 0.01
        if (Vsim not in dataV): # or Vsim==-30.0 or Vsim==-10.0 or Vsim==10.0 or Vsim==30.0:# or Vsim!=40.0:
            continue
        # else do something
        It = np.loadtxt(pCSim)
        data_idx = list(dataV).index(Vsim) + 1  # 0 is time
        #### match the time
        ItTW = ( (It[:,0]>2) & (It[:,0]<500) )
        timeData = [round(i*1e-3-50,3) for i in data[:,0]]
        ItTWData = [i for i, item in enumerate(It[ItTW,0]) if round(item,3) in set(timeData)]
        It = It[ItTW,:][ItTWData,:]
        
        dataTWIt = [i for i, item in enumerate(data[:,0]*1e-3-50) if round(item,3) in set(It[:,0])]
        data = data[dataTWIt,:]
        AllIt.append(It)
        AllData.append(data[:,data_idx])
    
    # define the objective function
    def f_iNa_g_tomin(x):
        (rto,rks,rkr,rncx,rnak,rna) = x
        ERR = 0.0
        for i in range(len(AllData)):
            I1 = AllIt[i][:,1]*rkr + AllIt[i][:,2]*rks + AllIt[i][:,3]*rto + AllIt[i][:,4]*rncx + AllIt[i][:,5]*rnak + AllIt[i][:,6]*rna
            I2 = AllData[i][:]
            ERR += np.linalg.norm(I1-I2)
        return ERR
    
    # the starting point
    x0 = [1., 10., 7., 5., 10., 1.]
    
    # the bounds
    xmin = [0., 0., 0., 0., 0., 0.] #[ 0., 18.55682001,  0. ,  0. , 0.78500075]#
    xmax = [100., 30., 20., 11., 50., 11.]
    
    # rewrite the bounds in the way required by L-BFGS-B
    bounds = [(low, high) for low, high in zip(xmin, xmax)]
    
    # use method L-BFGS-B because the problem is smooth and bounded
    minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds)
    res = basinhopping(f_iNa_g_tomin, x0, minimizer_kwargs=minimizer_kwargs)
    
    print "global minimum: x = ", res.x, ", f(x0) = ", res.fun
    (rtom,rksm,rkrm,rncxm,rnakm,rnam) = res.x
    # Reload data for plotting
    data = np.loadtxt(pathToData+dataFile)
    dataV = np.loadtxt(voltageFile)
    colors = get_spaced_colors(len(dataV))
    ccount = 0
    for pCSim in PCSims:
        Vsim = float(re.findall(pCSimReV,pCSim)[0]) - 0.01
        if (Vsim not in dataV) or Vsim==-30.0 or Vsim==-10.0 or Vsim==10.0 or Vsim==30.0:# or Vsim!=40.0:
            continue
        # else do something
        It = np.loadtxt(pCSim)
        data_idx = list(dataV).index(Vsim) + 1  # 0 is time
        I1m = It[:,1]*rkrm + It[:,2]*rksm + It[:,3]*rtom + It[:,4]*rncxm + It[:,5]*rnakm + It[:,6]*rnam
        
        
        # It[:,na_idx]+It[:,3]+It[:,1]*20+It[:,2]*10+It[:,4]+It[:,5];   It[:,3]*18;   It[:,1]*30
        ### It[:,3]*8+It[:,na_idx]+It[:,2]*10+It[:,1]*10   +   It[:,4]*4+It[:,5]*5
        #plt.plot(It[:,0],It[:,3]*7+It[:,na_idx]+It[:,2]*10   +   It[:,4]*2,color=colors[ccount],linewidth=2.0)#,label='V=%gmV'%Vsim)
        plt.plot(It[:,0], I1m, color=colors[ccount],linewidth=2.0)
        plt.plot(data[:,0]*1e-3-50, data[:,data_idx],color=colors[ccount])#,label='V=%gmV'%Vsim)
        ccount += 1
    plt.xlabel("t [ms]")
    plt.ylabel("I [A/F]")
    plt.xlim([-30, 530])
    plt.ylim([-10, 15])
    plt.title(dataFile)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3, fancybox=True, shadow=True)
    if save==None:
        plt.show()
    else:
        plt.savefig("%s_%s.png"%(save,dataFile[:-4]),bbox_inches='tight')
    plt.close('all')






def compareKTracesDiffAll(PCSims, dataFile, pathToData,save=None):
    ## Compare all voltage
    # load data for minimisation and format correction
    data = np.loadtxt(pathToData+dataFile)
    dataV = np.loadtxt(voltageFile)
    AllIt = []
    AllData = []
    for pCSim in PCSims:
        Vsim = float(re.findall(pCSimReV,pCSim)[0]) - 0.01
        if (Vsim not in dataV): # or Vsim==-30.0 or Vsim==-10.0 or Vsim==10.0 or Vsim==30.0:# or Vsim!=40.0:
            continue
        # else do something
        It = np.loadtxt(pCSim)
        data_idx = list(dataV).index(Vsim) + 1  # 0 is time
        #### match the time
        ItTW = ( (It[:,0]>2) & (It[:,0]<500) )
        timeData = [round(i*1e-3-50,3) for i in data[:,0]]
        ItTWData = [i for i, item in enumerate(It[ItTW,0]) if round(item,3) in set(timeData)]
        It = It[ItTW,:][ItTWData,:]
        
        dataTWIt = [i for i, item in enumerate(data[:,0]*1e-3-50) if round(item,3) in set(It[:,0])]
        data = data[dataTWIt,:]
        AllIt.append(It)
        AllData.append(data[:,data_idx])
    
    # define the objective function
    def f_iNa_g_tomin(x):
        (rto,rks,rkr,rncx,rnak,rna,rk1,rf,rcal,rcap) = x
        ERR = 0.0
        for i in range(len(AllData)):
            I1 = AllIt[i][:,1]*rkr + AllIt[i][:,2]*rks + AllIt[i][:,3]*rto + AllIt[i][:,4]*rncx + AllIt[i][:,5]*rnak + AllIt[i][:,6]*rna + AllIt[i][:,7]*rk1 + AllIt[i][:,8]*rf + AllIt[i][:,9]*rcal + AllIt[i][:,10]*rcap
            I2 = AllData[i][:]
            ERR += np.linalg.norm(I1-I2)
        return ERR
    
    # the starting point
    #x0 = [1., 10., 7., 5., 10., 1.,  1.,1.,1.,1.]
    x0 = np.random.random(10)*np.array([50., 50., 20., 11., 50., 11., 11.,11.,11.,11.])
    print 'init guess: ', x0
    # the bounds
    xmin = [0., 0., 0., 0., 0., 0., 0.,0.,0.,0.] #[ 0., 18.55682001,  0. ,  0. , 0.78500075]#
    xmax = [100., 100., 20., 11., 50., 11., 11.,11.,11.,11.]
    
    # rewrite the bounds in the way required by L-BFGS-B
    bounds = [(low, high) for low, high in zip(xmin, xmax)]
    
    # use method L-BFGS-B because the problem is smooth and bounded
    minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds)
    res = basinhopping(f_iNa_g_tomin, x0, minimizer_kwargs=minimizer_kwargs)
    
    print "global minimum: x = ", res.x, ", f(x0) = ", res.fun
    (rtom,rksm,rkrm,rncxm,rnakm,rnam,rk1m,rfm,rcalm,rcapm) = res.x
    # Reload data for plotting
    '''data = np.loadtxt(pathToData+dataFile)
    dataV = np.loadtxt(voltageFile)
    colors = get_spaced_colors(len(dataV))
    ccount = 0
    for pCSim in PCSims:
        Vsim = float(re.findall(pCSimReV,pCSim)[0]) - 0.01
        if (Vsim not in dataV) or Vsim==-30.0 or Vsim==-10.0 or Vsim==10.0 or Vsim==30.0:# or Vsim!=40.0:
            continue
        # else do something
        It = np.loadtxt(pCSim)
        data_idx = list(dataV).index(Vsim) + 1  # 0 is time
        I1m = It[:,1]*rkrm + It[:,2]*rksm + It[:,3]*rtom + It[:,4]*rncxm + It[:,5]*rnakm + It[:,6]*rnam + It[:,7]*rk1m + It[:,8]*rfm + It[:,9]*rcalm + It[:,10]*rcapm
        if Vsim == 50.0:
            idx_t500 = list(It[:,0]).index(490.0)
            #print "I_i/Itotal|v=50mV : ", [It[idx_t500,3]*rtom/I1m[idx_t500], It[idx_t500,2]*rksm/I1m[idx_t500], It[idx_t500,1]*rkrm/I1m[idx_t500], It[idx_t500,4]*rncxm/I1m[idx_t500], It[idx_t500,5]*rnakm/I1m[idx_t500], It[idx_t500,6]*rnam/I1m[idx_t500], It[idx_t500,7]*rk1m/I1m[idx_t500], It[idx_t500,8]*rfm/I1m[idx_t500], It[idx_t500,9]*rcalm/I1m[idx_t500], It[idx_t500,10]*rcapm/I1m[idx_t500]]
        
        
        # It[:,na_idx]+It[:,3]+It[:,1]*20+It[:,2]*10+It[:,4]+It[:,5];   It[:,3]*18;   It[:,1]*30
        ### It[:,3]*8+It[:,na_idx]+It[:,2]*10+It[:,1]*10   +   It[:,4]*4+It[:,5]*5
        #plt.plot(It[:,0],It[:,3]*7+It[:,na_idx]+It[:,2]*10   +   It[:,4]*2,color=colors[ccount],linewidth=2.0)#,label='V=%gmV'%Vsim)
        plt.plot(It[:,0], I1m, color=colors[ccount],linewidth=2.0)
        plt.plot(data[:,0]*1e-3-50, data[:,data_idx],color=colors[ccount])#,label='V=%gmV'%Vsim)
        ccount += 1
    plt.xlabel("t [ms]")
    plt.ylabel("I [A/F]")
    plt.xlim([-30, 530])
    plt.ylim([-10, 15])
    plt.title(dataFile)
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3, fancybox=True, shadow=True)
    if save==None:
        plt.show()
    else:
        plt.savefig("%s_%s.png"%(save,dataFile[:-4]),bbox_inches='tight')
    plt.close('all')'''





import sys
sys.path.append("/users/cholei/Downloads/cma-1.1.7")
import cma


def compareKTracesDiffAll_cma(PCSims, dataFile, pathToData,save=None):
    ## Compare all voltage
    # load data for minimisation and format correction
    data = np.loadtxt(pathToData+dataFile)
    dataV = np.loadtxt(voltageFile)
    AllIt = []
    AllData = []
    for pCSim in PCSims:
        Vsim = float(re.findall(pCSimReV,pCSim)[0]) - 0.01
        if (Vsim not in dataV): # or Vsim==-30.0 or Vsim==-10.0 or Vsim==10.0 or Vsim==30.0:# or Vsim!=40.0:
            continue
        # else do something
        It = np.loadtxt(pCSim)
        data_idx = list(dataV).index(Vsim) + 1  # 0 is time
        #### match the time
        ItTW = ( (It[:,0]>2) & (It[:,0]<500) )
        timeData = [round(i*1e-3-50,3) for i in data[:,0]]
        ItTWData = [i for i, item in enumerate(It[ItTW,0]) if round(item,3) in set(timeData)]
        It = It[ItTW,:][ItTWData,:]
        
        dataTWIt = [i for i, item in enumerate(data[:,0]*1e-3-50) if round(item,3) in set(It[:,0])]
        data = data[dataTWIt,:]
        AllIt.append(It)
        AllData.append(data[:,data_idx])
    
    # define the objective function
    def f_iNa_g_tomin(x):
        (rto,rks,rkr,rncx,rnak,rna,rk1,rf,rcal,rcap) = x
        ERR = 0.0
        for i in range(len(AllData)):
            I1 = AllIt[i][:,1]*rkr + AllIt[i][:,2]*rks + AllIt[i][:,3]*rto + AllIt[i][:,4]*rncx + AllIt[i][:,5]*rnak + AllIt[i][:,6]*rna + AllIt[i][:,7]*rk1 + AllIt[i][:,8]*rf + AllIt[i][:,9]*rcal + AllIt[i][:,10]*rcap
            I2 = AllData[i][:]
            ERR += np.linalg.norm(I1-I2)
        return ERR
    
    # the starting point
    #x0 = [1., 10., 7., 5., 10., 1.,  1.,1.,1.,1.]
    #x0 = np.random.random(10)*np.array([50., 50., 20., 50., 50., 11., 11.,11.,11.,11.])
    x0 = np.random.random(10)*70.
    print 'init guess: ', x0
    # the bounds
    xmin = [0., 0., 0., 0., 0., 0., 0.,0.,0.,0.] #[ 0., 18.55682001,  0. ,  0. , 0.78500075]#
    #xmax = [100., 100., 20., 11., 50., 11., 11.,11.,11.,11.]
    xmax = np.ones(10)*100.
    
    # cma.CMAOptions()  # returns all possible options
    options = { 'verb_plot': 1, 'seed':1234, 'verb_time':0, 'bounds': [xmin,xmax]}#'CMA_diagonal':100,
    #bounds = [(low, high) for low, high in zip(xmin, xmax)]
    
    # use method cma fmin
    res = cma.fmin(f_iNa_g_tomin, x0, 2.5, options)
    plt.show()
    #print "global minimum: x = ", res[0], ", f(x0) = ", res[1]
    #toCheck = [  0.  ,        51.81601703 ,  0.     ,      6.60383406  , 0.     ,      0.85861491 ,  0.       ,    0.     ,      2.20483769 ,  0.        ]
    #toCheck = [  3.73318154 , 19.68358847 ,  0.    ,       3.03944831 ,  0.     ,      1.22272337,   0.    ,       0.       ,    0.24469998  , 0.        ] 
    #toCheck = [ 0.     ,     7.61213826,  0.     ,     1.62430275  ,0.     ,     1.20509787  ,3.56607215 , 0.      ,    1.3819038  , 0.        ] 
    toCheck = [ 0.      ,    8.06918087,  2.70293071 , 5.10782478 , 0.  ,        0.85280746  ,1.85183471 , 0.    ,      0.    ,      0.        ] 
    for ii in range(10):
        if np.abs(res[0][ii]-toCheck[ii])>1e-3:
            raise Exception("Error!!! index %d does not match: (%g,%g)"%(ii,res[0][ii],toCheck[ii]))
        else:
            print "PASS!"
    (rtom,rksm,rkrm,rncxm,rnakm,rnam,rk1m,rfm,rcalm,rcapm) = res[0]
    # Reload data for plotting
    '''data = np.loadtxt(pathToData+dataFile)
    dataV = np.loadtxt(voltageFile)
    colors = get_spaced_colors(len(dataV))
    ccount = 0
    for pCSim in PCSims:
        Vsim = float(re.findall(pCSimReV,pCSim)[0]) - 0.01
        if (Vsim not in dataV) or Vsim==-30.0 or Vsim==-10.0 or Vsim==10.0 or Vsim==30.0:# or Vsim!=40.0:
            continue
        # else do something
        It = np.loadtxt(pCSim)
        data_idx = list(dataV).index(Vsim) + 1  # 0 is time
        I1m = It[:,1]*rkrm + It[:,2]*rksm + It[:,3]*rtom + It[:,4]*rncxm + It[:,5]*rnakm + It[:,6]*rnam + It[:,7]*rk1m + It[:,8]*rfm + It[:,9]*rcalm + It[:,10]*rcapm
        if Vsim == 50.0:
            idx_t500 = list(It[:,0]).index(490.0)
            #print "I_i/Itotal|v=50mV : ", [It[idx_t500,3]*rtom/I1m[idx_t500], It[idx_t500,2]*rksm/I1m[idx_t500], It[idx_t500,1]*rkrm/I1m[idx_t500], It[idx_t500,4]*rncxm/I1m[idx_t500], It[idx_t500,5]*rnakm/I1m[idx_t500], It[idx_t500,6]*rnam/I1m[idx_t500], It[idx_t500,7]*rk1m/I1m[idx_t500], It[idx_t500,8]*rfm/I1m[idx_t500], It[idx_t500,9]*rcalm/I1m[idx_t500], It[idx_t500,10]*rcapm/I1m[idx_t500]]
        
        
        # It[:,na_idx]+It[:,3]+It[:,1]*20+It[:,2]*10+It[:,4]+It[:,5];   It[:,3]*18;   It[:,1]*30
        ### It[:,3]*8+It[:,na_idx]+It[:,2]*10+It[:,1]*10   +   It[:,4]*4+It[:,5]*5
        #plt.plot(It[:,0],It[:,3]*7+It[:,na_idx]+It[:,2]*10   +   It[:,4]*2,color=colors[ccount],linewidth=2.0)#,label='V=%gmV'%Vsim)
        plt.plot(It[:,0], I1m, color=colors[ccount],linewidth=2.0)
        plt.plot(data[:,0]*1e-3-50, data[:,data_idx],color=colors[ccount])#,label='V=%gmV'%Vsim)
        ccount += 1
    plt.xlabel("t [ms]")
    plt.ylabel("I [A/F]")
    plt.xlim([-30, 530])
    plt.ylim([-10, 15])
    plt.title(dataFile)
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3, fancybox=True, shadow=True)
    if save==None:
        plt.show()
    else:
        plt.savefig("%s_%s.png"%(save,dataFile[:-4]),bbox_inches='tight')
    plt.close('all')'''








#print "Pure IKr and IKs; Imax for (data,IKr,IKs) = ", getMaxCurrent0(PCSims, dataFile, pathToData, fixVat50=False)
## Test the getMaxCurrent() is working
'''
Imax = []
for rr in np.arange(0,1.1,0.1):
    #print "PercenIKs=%g, Imax for data and sim: "%(rr*100), getMaxCurrent(PCSims, dataFile, pathToData, rr)
    Imax.append(getMaxCurrent(PCSims, dataFile, pathToData, rr, fixVat50=False)[1])

plt.plot(np.arange(0,1.1,0.1)*100, Imax, 'o-')
plt.xlabel("percentage of IKs")
plt.ylabel(r"$I_{max}$ [A/F]")
plt.title("0 percent = all IKr")
#plt.savefig("IKr_IKs_match_50mV.png")
plt.show()
'''

## Plot the current traces at a given ratio and setting
'''
for rr in np.arange(0,1.1,0.1):
    getCurrentTraces(PCSims, dataFile, pathToData, rr, fixVat50=False,save='rawTraces')
'''

## Plot iNa raw traces; data vs sim
"""
compareNaTraces(PCSims, dataFile, pathToData)
"""
'''
## Find composition of current in iNa raw traces
for i in range(1):
    compareKTracesDiffAll_cma(PCSims, dataFile, pathToData, None) #"globalmin")
plt.show()

'''


