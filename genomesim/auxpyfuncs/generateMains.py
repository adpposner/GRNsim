#!/home/russell/anaconda3/bin/python
from itertools import product
from configsetupnewnets import generateNetworkWorker,getMicroInt,doMain,cmdrun
from multiprocessing import Pool,cpu_count,Value
import os
nProc = 8#cpu_count()-2
pool = Pool(nProc)#,cores_idx=list(range(0,nProc)))

# def runSimsForNMicro(pth,nSims,nMirs):
#     pargs=[[pth,nSims,nMirs]]*nProc
#     [pargs[i].append(i) for i in range(nProc)]
#     pool.map(runWorker,pargs)

ctr= Value('i',0)

def runAllNetworksAtOnce(networkslist,simsPerProc):
    ival = getMicroInt()
    arglist=list(product(networkslist,ival,(3,),(2,)))
    pool.map(doMain,arglist)


def generateKineticExecutables(txckhlist=(1.0,),txcprodelist=(3.0,),txlkhlist=(10.0,),txlexplist=(3.0,),deckhlist=(6.0,),deckexplist=(3.0,),forcesdecay=(True,)):
    ptot=list(product(txckhlist,txcprodelist,txlkhlist,txlexplist,deckhlist,deckexplist,forcesdecay))
    totlen=len(ptot)
    for i in range(totlen):
        ptemp=list(ptot[i])
        ptemp.append(totlen)
        ptot[i]=tuple(ptemp)
    pool.map(generateNetworkWorker,ptot)
    
import numpy as np
if __name__ == "__main__":
    # scale = 100
    khtxc=(1.0,)#np.arange(1,7,5)
    khtxl=(10.0,)#np.arange(1,22,10)
    khdk=(6.0,)#np.arange(1,22,10)
    exr = (3.0,)#np.arange(1,4,1)
    # ptot=list(product(khtxc,exr,khtxl,exr,khdk,exr))
    # print(ptot)
    # nets=[]
    #generateKineticExecutables(khtxc,exr,khtxl,exr,khdk,exr)
    #generateKineticExecutables(khtxc,exr,khtxl,exr,khdk,exr)
    netlist=[]
    for r,d,f in os.walk('output'):
        rdr=os.path.split(r)
        if len(rdr) ==2:
            if "ngs" in rdr[1]:
                if "genes.json" in f:
                    netlist.append(r)
    print(netlist)
    runAllNetworksAtOnce(netlist,None)
