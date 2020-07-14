import sys
import numpy as np
import scipy.interpolate as interp

if __name__ == "__main__":
    workingDirectory   = sys.argv[1]
    (KTpts, KTwts)     = np.loadtxt(workingDirectory + '/KT.dat').T
    (KPHIpts, KPHIwts) = np.loadtxt(workingDirectory + '/KPHI.dat').T

    nMom     = 2
    nKT      = len(KTpts)
    nKPHI    = len(KPHIpts)

    dataFile = workingDirectory \
               + '/total_Pion_+_dN_dypTdpTdphi_allmoms_ev_no_df.dat'
    data     = np.loadtxt(dataFile, useCols=(2,4)).reshape([nKT,nKPHI,nMom])
    
    KT_integrated_moments      = np.dot(KTwts, KTpts * data)
    KPHI_integrated_moments    = np.tensordot(KPHIwts, data, (0,1))
    KT_KPHI_integrated_moments = np.dot(KPHIwts, KT_integrated_moments)
    
    moments_at_KT = interp.interp1D(KTpts, KPHI_integrated_moments, kind='linear')
    #moments_at_KT = interp.interp1D(KTpts, KPHI_integrated_moments, kind='quadratic')
    #moments_at_KT = interp.interp1D(KTpts, KPHI_integrated_moments, kind='cubic')
    
    KTvalue = 0.25 # GeV
    (x2, y2) = moments_at_KT(KTvalue)
    print(x2 + y2)