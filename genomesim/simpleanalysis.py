#!/usr/bin/python3
import numpy as np
import pandas as pd
import sys


#quick time-to-diff using fft approx of DHT - run with argv[1] as the prot file name
if __name__ == '__main__':
    q=pd.read_csv(sys.argv[1],sep='\t').values[:,1:]
    qq=q/np.linalg.norm(q,axis=0,keepdims=True) 
    a = np.sum(qq,axis=1)
    ft = np.fft.fft(a)
    lim = a.shape[0] //2 + 1
    ft[lim:] = 0
    ft[:lim] = (2.0 +0.j) * ft[:lim]
    ift=np.fft.ifft(ft)
    print("time to diff for {0} = {1}".format(sys.argv[1],np.argmin(np.imag(ift))))

