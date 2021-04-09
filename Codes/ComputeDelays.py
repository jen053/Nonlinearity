from io import StringIO
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, filtfilt, hilbert
import os.path
import time
from scipy.optimize import leastsq
from matplotlib import lines as mlines
import glob

class timearraydata:

    def __init__(self, name,nt,nch):
        self.name = name
        self.time = np.zeros(nt)    # creates a new empty list
        self.data = np.zeros([nch,nt])
        
class timedata:

    def __init__(self, name):
        self.name = name
        self.time = np.array    # creates a new empty list
        self.data = np.array
        
class shifts:

    def __init__(self, name):
        self.name = name
        self.raw = np.array    # creates a new empty list
        self.remove_noise = np.array    # creates a new empty list
        self.remove_noise_trend = np.array    # creates a new empty list
        self.remove_noise_wiggles = np.array    # creates a new empty list
        self.delays=np.array

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a
   
   
def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    #y = lfilter(b, a, data)
    y=filtfilt(b,a,data)
    return y

def butter_lowpass(highcut, fs, order=5):
    nyq = 0.5 * fs
    #print( "nyq is", nyq)
    high = highcut / nyq
    b, a = butter(order, [high])
    return b, a
   
   
def butter_lowpass_filter(data, highcut, fs, order=5):
    b, a = butter_lowpass(highcut, fs, order=order)
    #y = lfilter(b, a, data)
    y=filtfilt(b,a,data)
    return y


def readmultichannellabdat(fname,verbose=False):
# read in data from the scope and make a matrix of values from it

    if verbose is True:
        print(fname)
    tmp=np.loadtxt(fname,skiprows=2,delimiter=',')
    [nt,nch]=np.shape(tmp)
    nch=nch-1

    dataset=timearraydata('dataset',nt,nch)

    dataset.time=tmp[:,0]

    for ch in range(nch):
        dat1=np.transpose(tmp[:,ch+1])
        dataset.data[ch,:]=dat1
        
    return dataset

def readlabdat(fname,dataset,nch=1):
# read in data from the scope and make a matrix of values from it
    
    
    f = open(fname,'r')
    header1=f.readline()
    header2=f.readline()
    tmp1=[]
    tmp2=[]
    
    #

    for line in f:
        columns = line.split(',')
        if columns[1]=='\n':
            continue
        tmp1.append(float(columns[0]))
        tmp2.append(float(columns[1]))

    dataset.time=np.asarray(tmp1)
    dataset.data=np.asarray(tmp2)
    f.close()

def writelabdat(fname,dataset,nch=1):
# read in data from the scope and make a matrix of values from it
    
    
    f = open(fname,'w')
    header1='written from python\n'
    header2=time.strftime("%d/%m/%Y")+'\n'
    tmp1=[]
    tmp2=[]
    
    #

    f.write(header1)
    f.write(header2)

    nt=np.size(dataset.time)

    for it in range(nt):
        f.write(str(dataset.time[it])+','+str(dataset.data[it])+'\n')

    f.close()
        
   
def Compute_single_delay(PUMPprobe,probe,PUMP,ch,winlen, getCorrels=False):
#Compute the delay between the perturbed and un-perturbed probe signals
    ch_ind=0

    #first subtract the pump signal from the perturbed probe signal
    size_PUMPprobe=np.shape(PUMPprobe[ch_ind].data)[0]
    size_PUMP=np.shape(PUMP[ch_ind].data)[0]
    
    if size_PUMPprobe != size_PUMP:
        if size_PUMPprobe < size_PUMP:
            PUMP[ch_ind].data=np.delete(PUMP[ch_ind].data,range(size_PUMPprobe-1,size_PUMP-1))
        else:
            PUMPprobe[ch_ind].data=np.delete(PUMPprobe[ch_ind].data,range(size_PUMP-1,size_PUMPprobe-1))
        
    #Perturbed_probe=np.subtract(PUMP_only[ch_ind].data,PUMP_probe[ch_ind].data)
    Perturbed_probe=np.subtract(PUMPprobe[ch_ind].data,PUMP[ch_ind].data)
    #Perturbed_probe=PUMP_only[ch_ind].data
    #print np.subtract(PUMPprobe[ch_ind].data,PUMP[ch_ind].data)
    #remove means
#    Perturbed_probe_0mean=Perturbed_probe-np.mean(Perturbed_probe)
#    probe_only_0mean=probe_only[ch_ind].data-np.mean(probe_only[ch_ind].data)
    
    #filter data
    dt=PUMPprobe[ch_ind].time[1]-PUMPprobe[ch_ind].time[0]
    fs=1.0/dt
    probe_filter_cut=2.5e6
    filter_order=5    
    
    Perturbed_probe_filtered=butter_lowpass_filter(Perturbed_probe,probe_filter_cut,fs,filter_order)
    probe_filtered=butter_lowpass_filter(probe[ch_ind].data,probe_filter_cut,fs,filter_order)
    
    Perturbed_probe_filtered=Perturbed_probe_filtered-np.mean(Perturbed_probe_filtered)
    probe_filtered=probe_filtered-np.mean(probe_filtered)
    
#    Perturbed_probe_filtered=Perturbed_probe-np.mean(Perturbed_probe)
    #probe_only_filtered=probe_only[ch_ind].data-np.mean(probe_only[ch_ind].data)
#    Perturbed_probe_filtered=Perturbed_probe
#    Perturbed_probe_filtered=PUMPprobe[ch_ind].data#Perturbed_probe
#    probe_filtered=probe[ch_ind].data

 
    
    #plt.plot(Perturbed_probe_filtered)
    #plt.plot(probe_filtered,'r')
    #plt.show()
    
    probe_argmax=np.argmax(probe_filtered)
    time_probemax=probe[ch_ind].time[probe_argmax]
    
    tmin_to_correlate=time_probemax-winlen/2
    tmax_to_correlate=time_probemax+winlen/2

    index_of_tmin=np.argmin(np.abs(tmin_to_correlate-probe[ch_ind].time))
    index_of_tmax=np.argmin(np.abs(tmax_to_correlate-probe[ch_ind].time))


    tvec=np.arange(index_of_tmin,index_of_tmax)
    
    ntimes=np.size(tvec)
    tapsize=10
    tap=cosinetaper(ntimes,tapsize)


    Correlation_of_probes=np.correlate(probe_filtered[tvec]*tap,Perturbed_probe_filtered[tvec]*tap,"full")
    
    index_of_timeshift=np.argmax(Correlation_of_probes)
    #print index_of_timeshift
    nt=len(tvec)
    #print "nt is", nt
    dt=probe[ch_ind].time[1]-probe[ch_ind].time[0]
    #print "dt is", dt
    tvec_correlation=np.linspace(-nt*dt,nt*dt,2*nt-1)
    #print tvec_correlation
    #t_shift=tvec_correlation[index_of_timeshift]
    
    
    tvec_lags_to_fit=tvec_correlation[index_of_timeshift-2:index_of_timeshift+2]
    
    Correlation_to_fit=Correlation_of_probes[index_of_timeshift-2:index_of_timeshift+2]
    Fit_to_peak=np.polyfit(tvec_lags_to_fit,Correlation_to_fit,2)
    fit = np.poly1d(Fit_to_peak)
#     xx = np.linspace(min(tvec_lags_to_fit), max(tvec_lags_to_fit), 1000) 
#     plt.plot(tvec_lags_to_fit,Correlation_to_fit)
#     plt.plot(xx, fit(xx))
#     plt.show()
  
    MaxPeak_from_interpolation=-Fit_to_peak[1]/2/Fit_to_peak[0]
    t_shift=MaxPeak_from_interpolation
    
    if getCorrels == True:
        return t_shift, fit, Fit_to_peak, Correlation_to_fit, tvec_lags_to_fit
    else:
        return t_shift
    
def Compute_Filter_All_Delays(delayvec,delay_scale_fact,missing_data,delaystep,path,fname_part2,fnameend,winlen,recompute,flipdelays=0,dofilter=0, getCorrels=False):

    delayfname=path+'Delays.txt'


    if recompute==0:
        if os.path.isfile(delayfname):
            dat=np.genfromtxt(delayfname)
            tshifts=shifts('tshifts')
            tshifts.raw=dat[:,1]
            tshifts.remove_noise=dat[:,2]
            tshifts.remove_noise_trend=dat[:,3]
            tshifts.remove_noise_wiggles=dat[:,4]
            tshifts.delays=dat[:,0]
            return tshifts
    
#    print delayvec
    missing_data_index=(missing_data-delayvec[0])/delaystep
    delayvec=np.delete(delayvec,missing_data_index)
    #print delayvec

    activechannels=[1]
    activechannelnames=[]
    for k in activechannels:
        activechannelnames.append('ch'+str(k))
    
    probe_only=[timedata(ch) for ch in activechannelnames]
    PUMP_only=[timedata(ch) for ch in activechannelnames]
    PUMP_probe=[timedata(ch) for ch in activechannelnames]

    j=0
    tshift=np.empty(np.shape(delayvec)[0],'float')
    fits=[]
    params=[]
    correls=[]
    lags=[]
    for delay in delayvec:
        #print(delay)
        #fname_pump=path+'pu'+fname_part2+str(int(delay*delay_scale_fact))+fnameend
        #fname_probe=path+'pr'+fname_part2+str(int(delay*delay_scale_fact))+fnameend
        #fname_both=path+'pp'+fname_part2+str(int(delay*delay_scale_fact))+fnameend
        fname_pump=path+'pu'+fname_part2+str(delay*delay_scale_fact)+fnameend
        fname_probe=path+'pr'+fname_part2+str(delay*delay_scale_fact)+fnameend
        fname_both=path+'pp'+fname_part2+str(delay*delay_scale_fact)+fnameend
        i=0
        for ch in activechannels:
            readlabdat(fname_both,PUMP_probe[i])
            readlabdat(fname_probe,probe_only[i])
            readlabdat(fname_pump,PUMP_only[i])
            i=i+1
        #print(PUMP_probe,probe_only,PUMP_only)
        if getCorrels == True:
            tshift[j], fit, param, correl, lag=Compute_single_delay(PUMP_probe,probe_only,PUMP_only,1,winlen,getCorrels=True)
            fits.append(fit)
            params.append(param)
            correls.append(correl)
            lags.append(lag)
        else:
            tshift[j]=Compute_single_delay(PUMP_probe,probe_only,PUMP_only,1,winlen)
        j=j+1

    fs=1.0/(delaystep*1.0e-6)
    
    tshifts=shifts('tshifts')
    tshifts.raw=-tshift*1.0e9
#    tshifts.remove_noise=butter_lowpass_filter(tshifts.raw,150.e3,fs,2)
    if len(tshifts.raw)<15 or dofilter is 0:
        tshifts.remove_noise=tshifts.raw
        tshifts.remove_noise_trend=tshifts.raw
        tshifts.remove_noise_wiggles=tshifts.raw
    else:
        
        tshifts.remove_noise=butter_lowpass_filter(tshifts.raw,250.e3,fs,2)
        tshifts.remove_noise_trend=butter_bandpass_filter(tshifts.raw,50.e3,150.e3,fs,2)
        #tshifts.remove_noise_trend=butter_bandpass_filter(tshifts.raw,70.e3,120.e3,fs,2)
        tshifts.remove_noise_wiggles=butter_lowpass_filter(tshifts.raw,50.e3,fs,2)  #50


    if flipdelays is 1:
        tshifts.delays=delayvec*-1.0
    else:
        tshifts.delays=delayvec
       
    dat=np.empty([np.size(tshifts.raw),5])
    dat[:,0]=tshifts.delays
    dat[:,1]=tshifts.raw
    dat[:,2]=tshifts.remove_noise
    dat[:,3]=tshifts.remove_noise_trend
    dat[:,4]=tshifts.remove_noise_wiggles

    #np.savetxt(delayfname,dat)
    if getCorrels == True:                
        return tshifts, fits, params, correls, lags
    else:
        return tshifts
    
def Compute_Filter_All_Delays_JIM(delayvec,delay_scale_fact,missing_data,delaystep,path,fname_part2,fnameend,winlen):

    missing_data_index=(missing_data-delayvec[0])/delaystep
    delayvec=np.delete(delayvec,missing_data_index)

    activechannels=[1]
    activechannelnames=[]
    for k in activechannels:
        activechannelnames.append('ch'+str(k))
    
    probe_only=[timedata(ch) for ch in activechannelnames]
    PUMP_only=[timedata(ch) for ch in activechannelnames]
    PUMP_probe=[timedata(ch) for ch in activechannelnames]

    j=0
    tshift=np.empty(np.shape(delayvec)[0],'float')

    for delay in delayvec:
        fname_pump=path+fname_part2+'PUMPalone'+str(int(delay*delay_scale_fact))+fnameend
        fname_probe=path+fname_part2+'probealone'+str(int(delay*delay_scale_fact))+fnameend
        fname_both=path+fname_part2+'Pptogether'+str(int(delay*delay_scale_fact))+fnameend
        i=0
        for ch in activechannels:
            readlabdat(fname_both,PUMP_probe[i])
            readlabdat(fname_probe,probe_only[i])
            readlabdat(fname_pump,PUMP_only[i])
            i=i+1
        tshift[j]=Compute_single_delay(PUMP_probe,probe_only,PUMP_only,1,winlen)
        j=j+1

    fs=1/(delaystep*1e-6)

    tshifts=shifts('tshifts')
    tshifts.raw=-tshift*1e9
#    tshifts.remove_noise=butter_lowpass_filter(tshifts.raw,150e3,fs,2)
    tshifts.remove_noise=butter_lowpass_filter(tshifts.raw,250e3,fs,2)
    tshifts.remove_noise_trend=butter_bandpass_filter(tshifts.raw,70e3,120e3,fs,2)
    tshifts.remove_noise_wiggles=butter_lowpass_filter(tshifts.raw,50.e3,fs,2)
    tshifts.delays=delayvec
                            
    return tshifts
    

def Compute_Delay_2Sigs(sig1,sig2,winlen):
    
    #filter data
    dt=sig1.time[1]-sig1.time[0]
    fs=1.0/dt
    probe_filter_cut=2.5e6
    filter_order=5    

    Perturbed_probe_filtered=sig1.data
    probe_filtered=sig2.data


    
    #plt.plot(Perturbed_probe_filtered)
    #plt.plot(probe_filtered,'r')
    #plt.show()
    
    probe_argmax=np.argmax(probe_filtered)
    time_probemax=sig1.time[probe_argmax]
    
    tmin_to_correlate=time_probemax-winlen/2
    tmax_to_correlate=time_probemax+winlen/2

    index_of_tmin=np.argmin(np.abs(tmin_to_correlate-sig1.time))
    index_of_tmax=np.argmin(np.abs(tmax_to_correlate-sig1.time))

    tvec=np.arange(index_of_tmin,index_of_tmax)
    Correlation_of_probes=np.correlate(probe_filtered[tvec],Perturbed_probe_filtered[tvec],"full")
    index_of_timeshift=np.argmax(Correlation_of_probes)

    nt=len(tvec)

    tvec_correlation=np.linspace(-nt*dt,nt*dt,2*nt-1)

    #t_shift=tvec_correlation[index_of_timeshift]
    
    
    tvec_lags_to_fit=tvec_correlation[index_of_timeshift-2:index_of_timeshift+2]
    Correlation_to_fit=Correlation_of_probes[index_of_timeshift-2:index_of_timeshift+2]

    Fit_to_peak=np.polyfit(tvec_lags_to_fit,Correlation_to_fit,2)

    MaxPeak_from_interpolation=-Fit_to_peak[1]/2/Fit_to_peak[0]
    t_shift=MaxPeak_from_interpolation
    return t_shift

def cosinetaper(numpts,tapperc):
    ntap=int(np.floor(numpts*tapperc/100.0))
    if ntap%2 == 0:
        ntap=ntap+1

    tap=np.hanning(ntap)
    a1=int(round((ntap-1)/2))
    tap1=tap[0:a1]
    tap2=tap[a1:ntap]

    tap_out=np.ones(numpts)
    tap_out[0:a1]=tap1
    tap_out[numpts-a1-1:numpts]=tap2

    return tap_out

def GeneratePathlist(rootdir,excludelist,pathlist):
    
    #get list of subdirectories of current directory
    subdirlist=next(os.walk(rootdir))[1]

    #if there are no subdirectories, then pathlist is complete so return
    if len(subdirlist)==0:
        pathlist.append(rootdir+"/")
        return pathlist
    
    #remove any directories that are on the input exclude list from the list of subdirectories of the current directory
    removelist=[]
    for text in excludelist: 
        for dirname in subdirlist:
            if text in dirname:
                removelist.append(dirname)
        for dirname in removelist:
            subdirlist.remove(dirname)
            removelist=[]
            
    #now we have to go through all of the subdirectories of the current directory and find all of its subdirectories
    for dirname in subdirlist:
        pathlist=GeneratePathlist(rootdir+'/'+dirname,excludelist,pathlist)
        
    return pathlist

#this extracts the numbers from the end of a .csv filename
def getnum(a):
    nl=len(a)
    testval=a[nl-6:nl-4]
    if testval.isdigit()==True:
        return int(testval)
    else:
        testval=a[nl-5:nl-4]
        if testval.isdigit() is True:
            return int(testval)
        else:
            testval=a[nl-6:nl-4]
            if testval.isdigit() is True:
                return int(testval)
            else:
                print('problem in getnum')

def plotPUMP(pathbase='',pathlist=''):
#This plots a the PUMP data and returns the maximum values of the recorded data as a function of paths in the path list. It also plots the data for the entire pathlist, but has no plt.show() so the figure will show up later.  
#pathbase is the base path, and should be just a single string
#pathlist is a list of subdirectories of pathbase that should contain a set of .csv files with data.  The data should be pp#.csv pu#.csv and pr#.csv, where pp are data with both the pump and the probe on, pr with just the probe and pu with just the pump.  The # is the delay between the pump and probe.

    fig1=plt.figure()
    maxvalpump=[]
    for pathend in pathlist:
        fname_root=pathbase+pathend
        pumpdat=timedata('pumpdat')
        readlabdat(fname_root+"STrans_B4.csv",pumpdat)
        pathsplit=pathend.split('/')
        labeltxt=pathend#pathsplit[len(pathsplit)-2]
        linesym='-'
        #if pathsplit[0] == "Sample2":
#        if pathsplit[len(pathsplit)-5] == "May4":
#            linesym='--'
#        else:
#            linesym='-'

     
        plt0, =plt.plot(pumpdat.time*1e6,pumpdat.data,label=labeltxt,ms=2,linestyle=linesym)
        #plt0, =plt.plot(pumpdat.time*1e6,pumpdat.data,label=labeltxt)
        mycol=plt.getp(plt0,'color')
        plt.setp(plt0,'markerfacecolor',mycol)
        plt.setp(plt0,'markeredgecolor',mycol)
#    maxrange=np.where(np.logical_and(110<pumpdat.time*1e6,pumpdat.time*1e6<120))
        maxrange=np.where(np.logical_and(50<pumpdat.time*1e6,pumpdat.time*1e6<150))
        maxvalpump.append(np.max(pumpdat.data[maxrange]))

    
    #figname=pathbase+figname_end
    #plt.xlim([50,150])
    plt.xlabel('time ($\mu$s)')
    plt.ylabel('Amplitude (V)')
    plt.legend(loc='upper left')#bbox_to_anchor=(0, 1), loc='upper left', ncol=1)#plt.legend(loc=2,bbox_to_anchor=(0.5,0))
    #plt.savefig(figname,bbox_inches='tight')
    #plt.show()
    
    return maxvalpump




def PlotDelays(pathlist,maxvalpump,winlen=60e-6,savefig=0,pathbase='',compbase='',recompute=0,missing_data=np.array([]),delaystep=1.0,delay_scale_fact=1,doslopes=0,plotcolor='O'):
#This computes the delays and plots: (i) the delays for all of the data (ii) the delays as a function of input voltage (ii) the wiggle amplitude (which is actually the amplitude of a best-fit sinusoid) (iii) the delays as a function of the maximum recorded pump amplitude and (iv) the wiggle amplitude as a function of the maximum recorded pump amplitude.  It has limited capacity to fit some of these data to a line, but that's definitely not debugged.  Inputs are:
#pathlist: list of directory names that contain the .csv files (see description for plotPUMP
#maxvalpump: output from plotPUMP, which is the maximum value of the pump for each filename.  Should have the same dimensions of the pathlist
#winlen=60e-6: window over which the cross-correlation is computed
#savefig=0: =1 saves figures in the directory that the data are held in, as well as at levels up from there
#pathbase='': can be used so that the pathlist need not contain some root directory, so the code will look for files in pathbase+pathlist[ii]
#compbase='': used in the filenames of figures it outputs, so that you can say what two data sets are being compared in the output
#recompute=0: =1 will recompute the delays (e.g. to change the window length), and save them to disk (delays.txt) in the same directory as the data.  If this is set to 0 it will only compute the delays if it doesn't find the file delays.txt
#missing_data=np.array([]): if you put numbers in here it will skip those delays, useful if there are bad data points
#delaystep=1.0: step between subsequent delays
#delay_scale_fact=1: If you've named the files with some scaling factor so ppX.csv is for delay x/10 then set the scale factor to 10
#doslopes=0: =1 computes the slope of the maximum delay as a function of the input voltage or the maximum of the pump amplitude
#plotcolor='O': Sets how the color of the different lines in the plot are determined.  'V' uses the amplitude to set the color, otherwise it's orientaiton (I think) or random.  
#The labels in the legend are SxOyTzAw where x is the sample number, y is the orientation, z is the trial and w is the amplitude

    fig2=plt.figure(figsize=(20,8))
    fig3=plt.figure()
    fig4=plt.figure()
    fig5=plt.figure()
    fig6=plt.figure()
    fig7=plt.figure()
    
    
    slopevec=np.zeros(np.shape(maxvalpump))             

#missing_data=np.array([])#12,21,22,25,26,29,43])

    plotscale=1.0#maxval[0]/maxval
    doscal=0
    aval=0.0
    ii=0
    jj=0
    maxvalvec=np.zeros([7,1])
    xvalvec=np.zeros([7,1])
    for pathend in pathlist:
       
        if doscal==0:
            scal=1.0
        else:
            scal=plotscale[ii]
        ii=ii+1
        path=pathbase+pathend
        fname_part2=''
        fnameend='.csv'
        tshifts2=shifts('tshifts2')

        path_split=path.split('/')
        nlevels=len(path_split)
        delayflipstring=path_split[nlevels-4]
        if delayflipstring.count('PUMP') is 1:
            flipdelays=1
            print("delays have been flipped for"+path+pathend)

        else:
            flipdelays=0

        
        
        #figure out the delays that are recorded in the directory; note this doesn't check for missing files, just finds the range
        filelist=sorted(glob.glob(path+"pp*"),key=getnum)

#The below does it an old way, there's a better way (I think) below that

#
#        print filelist
#        nf=len(filelist)
#        if nf==0:
#            print('problem with the filelist'+str(path))
#        firstfile=filelist[0]
#        firstfilename=firstfile.split('/')
#        ndirs=len(firstfilename)
#        firstfilename2=firstfilename[ndirs-1]
#        delay1=firstfilename2[2:4]
#        if delay1.isdigit is False:
#            delay1=firstfilename2[2:3]
#            if delay1.isdigit is False:
#                delay1=firstfilename2[3:4]
#                if delay1.isdigit is False:
#                    print("problem getting delayvec",path)
#                    
#        lastfile=filelist[nf-1]
#        lastfilename=lastfile.split('/')
#        ndirs=len(lastfilename)
#        lastfilename2=lastfilename[ndirs-1]
#        delay2=lastfilename2[2:4]
#        if delay2.isdigit is False:
#            delay2=lastfilename2[2:3]
#            if delay2.isdigit is False:
#                delay2=lastfilename2[3:4]
#                if delay2.isdigit is False:
#                    print("problem getting delayvec",path)
#        delayvec=np.arange(float(delay1),float(delay2),delaystep)


        dvec=[]
        for fname in filelist:
            filesplit=fname.split('/')
            ndirs=len(filesplit)
            firstfilename2=filesplit[ndirs-1]
            ndig=len(firstfilename2)
            delay1=firstfilename2[2:ndig-4]
            dvec.append(float(delay1))

        delayvec=np.sort(dvec)

        
        #Compute the time delays
        tshifts2=Compute_Filter_All_Delays(delayvec,delay_scale_fact,missing_data,delaystep,path,fname_part2,fnameend,winlen,recompute,flipdelays=flipdelays)

        #A huge shift is usually a problem in the xcorr, typically because the signal strength on the probe is small, this excludes those data
        maxshift=np.max(np.abs(tshifts2.raw))
#        if maxshift >100.0:
#            print("shift is too big "+str(maxshift)+">100.0")
#            print("path is: "+path)
#            print("delay is: "+str(tshifts2.delays[np.argmax(np.abs(tshifts2.raw))]))
#            continue
            
        #split the path
        pathsplit=pathend.split('/')
        
        Snum=0
        Onum=1
        Tnum=0
        for word in pathsplit:
            nw=len(word)
            if nw==0:
                continue
            if word[0:4]=='Samp':
                Snum=int(word[nw-1])
            elif word[0:4]=='Tria':
                Tnum=word[nw-1]
            elif word[0:4]=='Orie':  #Note that this takes the first number it finds in the filename
                found=False
                for lett in word:
                    if lett.isdigit()==True:
                        Onum=int(lett)
                        found=True
                if not found:
                    print("Can't find orientation",path)
                        
            elif word[nw-1]=='V':
                Anum=float(word[0:nw-2])
        
        labeltxt='S'+str(Snum)+'O'+str(Onum)+'T'+str(Tnum)+'A'+str(Anum)
        plt.figure(fig2.number)
        plt1, =plt.plot(tshifts2.delays,tshifts2.remove_noise*scal)
        mycol=plt.getp(plt1,'color')

      
        if Onum == 2:
            plotsymbol='^'
        else:
            plotsymbol='o'
        if plotcolor=='V':
            strval=str(Anum)
            if strval=='5.0':
                mycol='k'
            elif strval=='5.5':
                mycol='c'
            elif strval=='6.0':
                mycol='r'
            elif strval=='6.5':
                mycol='m'
            elif strval=='7.0':
                mycol='b'
            elif strval=='7.5':
                mycol='g'
            elif strval=='8.0':
                mycol='y'
            elif strval=='8.5':
                mycol='b'
            elif strval=='9.0':
                mycol='c'
            elif strval=='9.5':
                mycol='r'
            elif strval=='10.':
                mycol='k'
            plt.setp(plt1,'color',mycol)
            
    
        mssize=10
        if Snum==1:
            mycol2='r'
            aval=-0.15
            if plotsymbol=='^':
                plotsymbol='*'
                mssize=14
            else:
                plotsymbol='s'
        elif Snum==2:
            mycol2='b'
            aval=-0.05
        elif Snum==3:
            mycol2='m'
            aval=0.05
            if plotsymbol=='^':
                plotsymbol='*'
                mssize=14
            else:
                plotsymbol='s'
        else:
            mycol2='c'
            aval=0.15
        
        
        plt.figure(fig2.number)
        plt.plot(tshifts2.delays,tshifts2.raw*scal,plotsymbol,color=mycol,label=labeltxt)
        plt.plot(tshifts2.delays,tshifts2.remove_noise_wiggles*scal,'--',color=mycol)
        #plt.figure(fig2.number)
        #plt.plot(tshifts2.delays,tshifts2.remove_noise_trend*scal,color=mycol,label=labeltxt)
    
        maxval=np.max(tshifts2.remove_noise_wiggles)
        wiggmaxval=np.max(tshifts2.remove_noise_trend)
    
        if doslopes==1:
            maxvalvec[jj]=maxval
            xvalvec[jj]=Anum
            if Anum==10.0:
                fit=np.polyfit(xvalvec,maxvalvec,1)
                slopevec[ii-1]=fit[0]
                jj=0
            plt.figure(fig7.number)
            snum=pathsplit[0][6]
            plt.plot(float(snum),slopevec[ii-1],plotsymbol,color=mycol2,ms=mssize)
    
        #guess_mean=0.0
        #guess_std = wiggmaxval
        #guess_phase = 0.0
        # Define the function to optimize, in this case, we want to minimize the difference
        # between the actual data and our "guessed" parameters
        #optimize_func = lambda x: x[0]*np.sin(tshifts2.delays*2.0*np.pi*90.0e-3+x[1]) + x[2] - tshifts2.remove_noise_trend
        #est_std, est_phase, est_mean = leastsq(optimize_func, [guess_std, guess_phase, guess_mean])[0]

        #wiggmaxval=np.abs(est_std)
        #plt.figure()
        #plt.plot(tshifts2.delays,tshifts2.remove_noise_trend)
        #fitdat=est_std*np.sin(tshifts2.delays*2.0*np.pi*90.0e-3+est_phase) + est_mean
        #plt.plot(tshifts2.delays,fitdat)

        
        xval=Anum
        xval=xval+aval
    
        plt.figure(fig3.number)
        plt.plot(xval,maxval,plotsymbol,color=mycol2,ms=mssize)
        #plt.xlim([4,11])
        #plt.ylim([0,14])
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
    
        plt.figure(fig4.number)
        plt.plot(xval,wiggmaxval,plotsymbol,color=mycol2,ms=mssize)
        plt.ylim([0,4])
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
    
        plt.figure(fig5.number)
        plt.plot(maxvalpump[ii-1],maxval,plotsymbol,color=mycol2,ms=mssize)
        #plt.xlim([0,100])
        #plt.ylim([0,14])
      

        plt.figure(fig6.number)
        plt.plot(maxvalpump[ii-1],wiggmaxval,plotsymbol,color=mycol2,ms=mssize)
        plt.ylim([0,4])
 
    marker1 = mlines.Line2D([], [], linestyle='none',color='r', marker='*',ms=14, label='Sample 1, Orientation 4')
    marker2 = mlines.Line2D([], [], linestyle='none',color='r', marker='s',ms=10, label='Sample 1, Orientation 3')
    marker3 = mlines.Line2D([], [], linestyle='none',color='b', marker='^',ms=10, label='Sample 2, Orientation 2')
    marker4 = mlines.Line2D([], [], linestyle='none',color='b', marker='o',ms=10, label='Sample 2, Orientation 1')
    marker5 = mlines.Line2D([], [], linestyle='none',color='m', marker='*',ms=14, label='Sample 3, Orientation 4')
    marker6 = mlines.Line2D([], [], linestyle='none',color='m', marker='s',ms=10, label='Sample 3, Orientation 3')
    marker7 = mlines.Line2D([], [], linestyle='none',color='c', marker='^',ms=10, label='Sample 4, Orientation 2')
    marker8 = mlines.Line2D([], [], linestyle='none',color='c', marker='o',ms=10, label='Sample 4, Orientation 1')

    plt.figure(fig3.number)
    plt.legend(handles=[marker1,marker2,marker3,marker4,marker5,marker6,marker7,marker8],loc=2)
    plt.xlim([4,11])
        #plt.legend(loc=2)
    plt.xlabel('Pump Voltage (V)',fontsize=16)
    plt.ylabel('Max time delay (ns)',fontsize=16)
    plt.title('Comparison over all data',fontsize=16)
    plt.legend()

    plt.figure(fig4.number)
    #plt.legend(handles=[marker1,marker2,marker3,marker4,marker5,marker6,marker7,marker8],loc=1)
    plt.xlim([4,11])
    plt.xlabel('Pump Voltage (V)',fontsize=16)
    plt.ylabel('Wiggle Amplitude (ns)',fontsize=16)
    plt.title('Comparison over all data',fontsize=16)
    #plt.legend()

    plt.figure(fig5.number)
    #plt.legend(handles=[marker1,marker2,marker3,marker4,marker5,marker6,marker7,marker8],loc=2)
    plt.xlabel('Max Pump Amplitude (V)',fontsize=16)
    plt.ylabel('Max time delay (ns)',fontsize=16)
    plt.title('Comparison over all data',fontsize=16)
    #plt.legend()

    plt.figure(fig6.number)
    plt.legend(handles=[marker1,marker2,marker3,marker4,marker5,marker6,marker7,marker8],loc=1)
    plt.xlabel('Max Pump Amplitude (V)',fontsize=16)
    plt.ylabel('Wiggle Amplitude (ns)',fontsize=16)
    plt.title('Comparison over all data',fontsize=16)
    plt.legend()

    plt.figure(fig2.number)
    plt.legend(loc=2)
    plt.xlabel('phase delay (us)',fontsize=18)
    plt.ylabel('time delay (ns)',fontsize=18)
    plt.title('Raw data',fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    
    plt.legend()
    

    if savefig==1:
        
        figname_data=pathbase+'Compare_'+compbase+'NonLinDat.png'
        figname1=pathbase+'Compare_'+compbase+'.png'
        figname2=pathbase+'Compare_'+compbase+'_Wiggles.png'
        figname3=pathbase+'Compare_'+compbase+'_PumpAmp.png'
        figname4=pathbase+'Compare_'+compbase+'_PumpAmp_Wiggles.png'
        plt.figure(fig2.number)
        plt.savefig(figname_data,bbox_inches='tight')
 
        plt.figure(fig3.number)
        plt.savefig(figname1,bbox_inches='tight')
        plt.figure(fig4.number)
        plt.savefig(figname2,bbox_inches='tight')
        plt.figure(fig5.number)
        plt.savefig(figname3,bbox_inches='tight')
        plt.figure(fig6.number)
        plt.savefig(figname4,bbox_inches='tight')
    plt.show()


def TimeShift_TwoSigs(sig1,sig2,window,tvec):
    
    nt=np.size(window)
    dt=tvec[1]-tvec[0]
    tapsize=10

    tvec_correlation=np.linspace(-nt*dt,nt*dt,2*nt-1)

    tap=cosinetaper(nt,tapsize)

    
    sig1_tocorr=sig1[window]*tap
    sig2_tocorr=sig2[window]*tap
    
    plt.plot(tvec[window],sig1_tocorr)
    plt.plot(tvec[window],sig2_tocorr)

    Correlation_of_sigs=np.correlate(sig1_tocorr,sig2_tocorr,"full")
    
    index_of_timeshift=np.argmax(Correlation_of_sigs)
     
    tvec_lags_to_fit=tvec_correlation[index_of_timeshift-2:index_of_timeshift+2]
    Correlation_to_fit=Correlation_of_sigs[index_of_timeshift-2:index_of_timeshift+2]

    Fit_to_peak=np.polyfit(tvec_lags_to_fit,Correlation_to_fit,2)
    
   
    MaxPeak_from_interpolation=-Fit_to_peak[1]/2/Fit_to_peak[0]
    t_shift=MaxPeak_from_interpolation

    return t_shift

def GetDelays(pathlist,maxvalpump,winlen=60e-6,savefig=0,pathbase='',compbase='',recompute=0,missing_data=np.array([]),delaystep=1.0,delay_scale_fact=1,doslopes=0,plotcolor='O'):
#This computes the delays and plots: (i) the delays for all of the data (ii) the delays as a function of input voltage (ii) the wiggle amplitude (which is actually the amplitude of a best-fit sinusoid) (iii) the delays as a function of the maximum recorded pump amplitude and (iv) the wiggle amplitude as a function of the maximum recorded pump amplitude.  
#It has limited capacity to fit some of these data to a line, but that's definitely not debugged.  Inputs are:
#pathlist: list of directory names that contain the .csv files (see description for plotPUMP
#maxvalpump: output from plotPUMP, which is the maximum value of the pump for each filename.  Should have the same dimensions of the pathlist
#winlen=60e-6: window over which the cross-correlation is computed
#savefig=0: =1 saves figures in the directory that the data are held in, as well as at levels up from there
#pathbase='': can be used so that the pathlist need not contain some root directory, so the code will look for files in pathbase+pathlist[ii]
#compbase='': used in the filenames of figures it outputs, so that you can say what two data sets are being compared in the output
#recompute=0: =1 will recompute the delays (e.g. to change the window length), and save them to disk (delays.txt) in the same directory as the data.  If this is set to 0 it will only compute the delays if it doesn't find the file delays.txt
#missing_data=np.array([]): if you put numbers in here it will skip those delays, useful if there are bad data points
#delaystep=1.0: step between subsequent delays
#delay_scale_fact=1: If you've named the files with some scaling factor so ppX.csv is for delay x/10 then set the scale factor to 10
#doslopes=0: =1 computes the slope of the maximum delay as a function of the input voltage or the maximum of the pump amplitude
#plotcolor='O': Sets how the color of the different lines in the plot are determined.  'V' uses the amplitude to set the color, otherwise it's orientaiton (I think) or random.  
#The labels in the legend are SxOyTzAw where x is the sample number, y is the orientation, z is the trial and w is the amplitude

   
    shiftslist=[]
    
#missing_data=np.array([])#12,21,22,25,26,29,43])

    ii=0
    for pathend in pathlist:
       
        ii=ii+1
        path=pathbase+pathend
        fname_part2=''
        fnameend='.csv'
        tshifts2=shifts('tshifts2')

        path_split=path.split('/')
        nlevels=len(path_split)
        delayflipstring=path_split[nlevels-4]
        if delayflipstring.count('PUMP') is 1:
            flipdelays=1
            print("delays have been flipped for"+path+pathend)

        else:
            flipdelays=0

        
        
        #figure out the delays that are recorded in the directory; note this doesn't check for missing files, just finds the range
        filelist=sorted(glob.glob(path+"pp*"),key=getnum)



        dvec=[]
        for fname in filelist:
            filesplit=fname.split('/')
            ndirs=len(filesplit)
            firstfilename2=filesplit[ndirs-1]
            ndig=len(firstfilename2)
            delay1=firstfilename2[2:ndig-4]
            dvec.append(float(delay1))

        delayvec=np.sort(dvec)

        
        #Compute the time delays
        tshifts2=Compute_Filter_All_Delays(delayvec,delay_scale_fact,missing_data,delaystep,path,fname_part2,fnameend,winlen,recompute,flipdelays=flipdelays)
        shiftslist.append(tshifts2)
        
        
        

    return shiftslist
    
