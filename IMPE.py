from itertools import permutations
from math import log

def pe(TS, m, tstep=1):
    
    #List of all possible permutations
    m_list=[x for x in range(1, m+1)]
    l_perm=list(permutations(m_list))
    
    #Create dictionary of permutations (key) with number of matches (values) initialised at 0
    max_perm=len(l_perm) #total number of permutations
    m_perm=dict(zip(l_perm, [0 for x in range(1, max_perm+1)]))
    
    #loop across each value of the timeseries to obtain the d-sized vectors
    for t in range(1, len(TS)+1-tstep*(m-1)):
        index_start=t-1
        index_end=t-1+m
        
        v_sort=TS[index_start:index_end]
        v_sort.sort()
        
        v_perm=dict(zip(v_sort, m_list))
        
        v=TS[index_start:index_end]
        v_order=[v_perm[value] for value in v]
        
        #Look up and update dictionary with matches
        m_perm[tuple(v_order)]=m_perm[tuple(v_order)]+1

    #Calculate permutation entropy
    total=sum(m_perm.values())
    
    h_perm={x: -1* m_perm[x]/total * log(m_perm[x]/total) for x in m_perm if m_perm[x]!=0}
    H=sum(h_perm.values())
    
    return H


import numpy as np
def CoarseGrain(iSig, s):
    oSig=[]
    N=len(iSig)
    for i in list(range(1, int(N/s)+1)):
        oSig=oSig+[np.mean(iSig[(i-1)*s:i*s])]
    return oSig        


import numpy as np
def impe(TS, m, tstep, scale):
    
    if scale==1: 
        #calculate and return pe of the original TS, no coarsegraining required
        imPE=pe(TS, m, tstep)
        
        return imPE
    
    #Create temporary list to store pe values of s coarsegrain TSs
    PE_temp=[np.nan for x in range(1, scale+1)]
    
    #coarsegrain the signal, s coarsegrained TSs, calculate pe of each
    for i in list(range(0, scale)):
        TS_coarse=CoarseGrain(TS[i:], scale)
        
        #calculate pe for scale s for each coarsegrained TS
        H=pe(TS_coarse, m, tstep)
        
        #save value in temporary list
        PE_temp[i]=H
        
    #calculate impe for scale s, by taking mean of all pe values of each coarsegrained TS
    imPE=np.mean(PE_temp) #produces nan if permutation entropy cannot be obtained
    imPE_norm=imPE/log2(factorial(m))
    
    return imPE, imPE_norm