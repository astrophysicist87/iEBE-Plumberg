import sys
import numpy as np
import scipy.interpolate as interp
import matplotlib
import matplotlib.pyplot as plt

if __name__ == "__main__":
    workingDirectory   = sys.argv[1]
    (KTpts, KTwts)     = np.loadtxt(workingDirectory + '/KT.dat').T
    (KPHIpts, KPHIwts) = np.loadtxt(workingDirectory + '/KPHI.dat').T

    nMom     = 3
    nKT      = len(KTpts)
    nKPHI    = len(KPHIpts)

    dataFile = workingDirectory \
               + '/total_Pion_+_dN_dypTdpTdphi_allmoms_ev_no_df.dat'
    data     = np.loadtxt(dataFile).reshape([15,nKPHI,nKT])[[0,2,4]]
    
    KTdata = np.tile(KTpts,(nMom, nKPHI,1))
    
    KT_integrated_moments      = np.tensordot(KTwts, KTdata, (0,-1))
    KPHI_integrated_moments    = np.tensordot(KPHIwts, data, (0,1)) / (2.0*np.pi)
    KT_KPHI_integrated_moments = np.tensordot(KPHIwts, KT_integrated_moments, (0,-1)) / (2.0*np.pi)
    
    #moments_at_KT = interp.interp1d(KTpts, KPHI_integrated_moments, kind='linear')
    #moments_at_KT = interp.interp1d(KTpts, KPHI_integrated_moments, kind='quadratic')
    moments_at_KT = interp.interp1d(KTpts, KPHI_integrated_moments, kind='cubic')
    
    #print KPHI_integrated_moments
    
    KTvalue = 0.25 # GeV
    (norm, x2, y2) = moments_at_KT(KTvalue)
    print np.sqrt(0.5*(x2+y2)/norm)
    #print moments_at_KT(np.arange(0.01, 1.01, 0.1))
    
    
    '''fig, ax = plt.subplots()
    ax.scatter(KTpts, KPHI_integrated_moments[0],color='blue')
    ax.plot(np.arange(0.01, 1.01, 0.01), moments_at_KT(np.arange(0.01, 1.01, 0.01))[0],'b-')
    
    ax.scatter(KTpts, KPHI_integrated_moments[1],color='red')
    ax.plot(np.arange(0.01, 1.01, 0.01), moments_at_KT(np.arange(0.01, 1.01, 0.01))[1],'r-')
    
    ax.scatter(KTpts, KPHI_integrated_moments[2],color='green')
    ax.plot(np.arange(0.01, 1.01, 0.01), moments_at_KT(np.arange(0.01, 1.01, 0.01))[2],'g-')
    
    ax.set_xlim(0.0, 1.0)
    
    plt.show()'''