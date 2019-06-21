#!/home/russell/anaconda3/bin/python
import shutil
import random as rnd
from scipy.stats import rv_continuous
import xml.etree.ElementTree as et
import subprocess
import os
import timeit
from itertools import chain,repeat
from enum import Enum
from multiprocessing import Pool,cpu_count,Value,Lock
import numpy as np
from numbers import Number
from collections import Iterable


from collections import namedtuple
tfeffectxmlsearchname='.//TF/effect/continuous_distribution/distribution'
TFEffectsLogNorm=namedtuple('TFEffectsLogNorm',["mu","sigma","scale","displacement","type","weight"])
KineticsVars = namedtuple("KineticsVars",["txcK","txcExp","txlK","txlExp","decK","decExp","forcesdecay"])
xmlsearchnames=[".//TXC/K_ONE_HALF",".//TXC/PROD_EXPONENT",".//TXL/K_ONE_HALF",".//TXL/PROD_EXPONENT",
        ".//MESS_DECAY/K_ONE_HALF",".//MESS_DECAY/EXPONENT",".//kinetics/MESSENGERDECAYFORCESMICRODECAY"]
ConnectionDistribution=namedtuple('ConnectionDistribution',
["maxVal","minVal","successP","type","weight"])

class ContentType(Enum):
    FLOAT = 0
    INT = 1
    CONTINUOUSDISTRIBUTION = 2
    DISCRETEDISTRIBUTION = 3
    STRING = 4
    BOOLEAN = 7
    NONE = 8

def retNone(*args,**kwargs):
    return None
def boolFromString(s):
    if s.lower() == "true":
        return True
    else:
        return False


def clsInstanceFor(CTEN):
    if CTEN == ContentType.FLOAT:
        return float
    elif CTEN == ContentType.INT:
        return int
    elif CTEN == ContentType.CONTINUOUSDISTRIBUTION:
        return AbstractContinuousDistribution
    elif CTEN == ContentType.DISCRETEDISTRIBUTION:
        return AbstractDiscreteDistribution
    elif CTEN == ContentType.STRING:
        return str
    elif CTEN == ContentType.BOOLEAN:
        return bool
    # elif CTEN == ContentType.NONE:
    #     return None
    else:
        print(CTEN)
        raise NotImplementedError(CTEN.name,CTEN.value)

def clsCtorFor(CTEN):
    if CTEN == ContentType.FLOAT:
        return float
    elif CTEN == ContentType.INT:
        return int
    elif CTEN == ContentType.CONTINUOUSDISTRIBUTION:
        return ContinuousDistributionXMLFactory
    elif CTEN == ContentType.DISCRETEDISTRIBUTION:
        return DiscreteDistributionXMLFactory
    elif CTEN == ContentType.STRING:
        return str
    elif CTEN == ContentType.BOOLEAN:
        return boolFromString
    elif CTEN == ContentType.NONE:
        return retNone
    else:
        print(CTEN)
        raise NotImplementedError(CTEN.name,CTEN.value)


class LogNormTester:
    class LGNRMDIST(rv_continuous):
        def _pdf(self,x,mu,sigma,displacement,scal):
            if x>displacement:
                return np.exp(-((np.log((x-displacement)/scal)-mu)**2)/(2*(sigma**2)))/(sigma*(x-displacement)*np.sqrt(2*np.pi))
            else:
                return 0
    
    def __init__(self,mu,sigma,displacement,scal):
        pass

class AbstractConfigSection:
    __slots__= ()
    xmlSearchString = None

    #__SlotClasses__ = ()
    def __init__(self,**kwargs):
        for attribute in self.__slots__:
            setattr(self,attribute,None)
        if kwargs:
            self.fromDict(**kwargs)

    def __setattr__(self,key,value):
        ind = self.__slots__.index(key)
        slttype = self.__SlotClasses__[ind]
        if value is not None:
            #print(key,value)
            kls = clsCtorFor(slttype)
            if kls == int or kls == float:
                #print(kls)
                if not isinstance(value,Number):
                    raise AttributeError("value " + str(type(value)) + " for key " + str(key) +
                    " not a number" )
            elif slttype == ContentType.BOOLEAN:
                if not isinstance(value,bool):
                    raise AttributeError("value " + str(type(value)) + " for key " + str(key) +
                    " not a boolean" )
            elif not isinstance(value,clsInstanceFor(slttype)):
                raise AttributeError("value " + str(type(value)) + " for key " + str(key) +
                " not a " + str(kls))
            
        super().__setattr__(key,value)
    def __getitem__(self, key):
        return getattr(self, key)
    
    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __len__(self):
        return len(self.__slots__)

    def __iter__(self):
        return iter(self.__slots__)

    def items(self):
        for attribute in self.__slots__:
            yield attribute, getattr(self, attribute)
    
    def fromDict(self,**kwargs):
        for k,v in kwargs.items():
            setattr(self,k,v)

    @classmethod
    def getSection(kls,et):
        return et.find(kls.xmlSearchString)

    @classmethod
    def generateSweep(cls,**kwargs):
        prods=[]
        for key,value in kwargs.items():
            if key not in cls.__slots__:
                raise AttributeError("key " + str(key) + " not in slots for " + str(cls) + str(cls.__slots__))
            if not isinstance(value,Iterable):
                #print(key,value)
                raise AttributeError(str(value) +" for " + str(key) + " is not iterable!")
            prods.append(list(zip(repeat(key,len(value)),value)))
        for item in product(*prods):
            yield item
  
    def getTextValue(self,sec,str):
        return sec.find(str).text

    def getNumeric(self,sec,str):
        return float(sec.find(str))

    @classmethod
    def fromETree(cls,et):
        raise NotImplementedError(str(cls))

    def toETree(self,et):
        raise NotImplementedError(str(self.__class__))

class AbstractContinuousDistribution(AbstractConfigSection):
    __slots__ = ("mu","scale","sigma","displacement","type","weight")
    __SlotClasses__ = (ContentType.FLOAT,ContentType.FLOAT,ContentType.FLOAT,ContentType.FLOAT,ContentType.STRING,ContentType.FLOAT)
    xmlSearchString = "continuous_distribution/distribution"

    @classmethod
    def fromETree(cls,et):
        loc=cls.getSection(et)
        d={}
        for item,value in loc.attrib.items():
            try:
                if item not in cls.__slots__:
                    raise AttributeError("item " + str(item) + " not in slots for Continuous Distn")
                kls = clsCtorFor(cls.__SlotClasses__[cls.__slots__.index(item)])
                d[item]=kls(value)
            except Exception as e:
                print(item,value)
                raise e
        return cls(**d)

    def toETree(self,et):
        loc=self.__class__.getSection(et)
        for fieldname,value in self.items():
            if isinstance(value,str):
                loc.set(fieldname,value)
            elif isinstance(value,Number):
                loc.set(fieldname,str(value))
            else:
                raise ValueError(fieldname+ " " + value)

class UniformContinuousDistribution(AbstractContinuousDistribution):
    __slots__ = ("minVal","maxVal","type","weight")
    __SlotClasses__ = (ContentType.FLOAT,ContentType.FLOAT,ContentType.STRING,ContentType.FLOAT)

class LogNormalContinuousDistribution(AbstractContinuousDistribution):
    __slots__ = ("mu","scale","sigma","displacement","type","weight")
    __SlotClasses__ = (ContentType.FLOAT,ContentType.FLOAT,ContentType.FLOAT,ContentType.FLOAT,ContentType.STRING,ContentType.FLOAT)



def ContinuousDistributionXMLFactory(et):
    loc = AbstractContinuousDistribution.getSection(et)
    if "type" not in loc.attrib:
        raise AttributeError("Type of continuous distn not specified " + str(loc))
    disttype = loc.attrib["type"].lower()
    if disttype == "uniform":
        return UniformContinuousDistribution.fromETree(et)
    elif disttype == "lognormal":
        return LogNormalContinuousDistribution.fromETree(et)
    else:
        raise ValueError("No continuous distribution of type "+ disttype)


class AbstractDiscreteDistribution(AbstractConfigSection):
    # argsdict = {
    #     "binomial": {"minVal","maxVal","successP","nSuccess","weight"},
    #     "negbinomial": {"minVal,maxVal","successP","nSuccess",}
    # }
    __slots__ = ("minVal","maxVal","successP","nSuccess","type","weight")
    __SlotClasses__ = (ContentType.INT,ContentType.INT,ContentType.FLOAT,ContentType.FLOAT,ContentType.STRING,ContentType.FLOAT)
    xmlSearchString = "discrete_distribution/distribution"

    @classmethod
    def fromETree(cls,et):
        loc=cls.getSection(et)
        d={}
        for item,value in loc.attrib.items():
            if item not in cls.__slots__:
                raise AttributeError("item " + str(item) + " not in slots for " + str(cls))
            kls = clsCtorFor(cls.__SlotClasses__[cls.__slots__.index(item)])
            d[item]=kls(value)
        return cls(**d)

    def toETree(self,et):
        loc=self.__class__.getSection(et)
        for fieldname,value in self.items():
            if isinstance(value,str):
                loc.set(fieldname,value)
            elif isinstance(value,Number):
                loc.set(fieldname,str(value))
            else:
                raise ValueError(fieldname+ " " + value)

class BinomialDiscreteDistribution(AbstractDiscreteDistribution):
    __slots__ = ("minVal","maxVal","successP","type","weight")
    __SlotClasses__ = (ContentType.INT,ContentType.INT,ContentType.FLOAT,ContentType.STRING,ContentType.FLOAT)

class NegativeBinomialDiscreteDistribution(AbstractDiscreteDistribution):
    __slots__ = ("minVal","maxVal","successP","nSuccess","type","weight")
    __SlotClasses__ = (ContentType.INT,ContentType.INT,ContentType.FLOAT,ContentType.FLOAT,ContentType.STRING,ContentType.FLOAT)

    

def DiscreteDistributionXMLFactory(et):
    loc = AbstractDiscreteDistribution.getSection(et)
    if "type" not in loc.attrib:
        raise AttributeError("Type of discrete distn not specified " + str(loc))
    disttype = loc.attrib["type"].lower()
    if disttype == "binomial":
        return BinomialDiscreteDistribution.fromETree(et)
    elif disttype == "negbinomial":
        return NegativeBinomialDiscreteDistribution.fromETree(et)
    else:
        raise ValueError("No discrete distribution of type "+ disttype)

class TFBaseRates(AbstractConfigSection):
    __slots__ = ("production","decay","affinity","frequency","effect")
    __SlotClasses__=(ContentType.FLOAT,ContentType.FLOAT,ContentType.FLOAT,
    ContentType.FLOAT,ContentType.CONTINUOUSDISTRIBUTION)
    xmlSearchString="baseRates/TF"

    @classmethod
    def fromETree(cls,et):
        loc = cls.getSection(et)
        d={}
        for item,tp in zip(cls.__slots__,cls.__SlotClasses__):
            kls = clsCtorFor(tp)
            #print(kls,tp)
            l=loc.find(item)
            #print("hello",kls,l.text)
            if issubclass(clsInstanceFor(tp),AbstractContinuousDistribution):
                d[item]=kls(l)
            elif issubclass(clsInstanceFor(tp),AbstractConfigSection):
                c[item]=kls.fromETree(l)
            else:
                d[item] = kls(l.text)
        #print(d)
        return cls(**d)

    def toETree(self,et):
        loc=self.__class__.getSection(et)
        for fieldname,value in self.items():
            if isinstance(value,str):
                loc.find(fieldname).text=value
            elif isinstance(value,Number):
                loc.find(fieldname).text=str(value)
            elif isinstance(value,AbstractConfigSection):
                value.toETree(loc.find(fieldname))
            else:
                raise ValueError(str(fieldname) + " " +str(value))



class MessBaseRates(AbstractConfigSection):
    __slots__ = ("production","decay")
    __SlotClasses__=(ContentType.FLOAT,ContentType.FLOAT)
    xmlSearchString="baseRates/mess"

    @classmethod
    def fromETree(cls,et):
        loc = cls.getSection(et)
        d={}
        for item,tp in zip(cls.__slots__,cls.__SlotClasses__):
            kls = clsCtorFor(tp)
            l=loc.find(item)
            #print("hello",kls,l.text)
            if issubclass(kls,AbstractConfigSection):
                d[item]=kls.fromETree(l)
            else:
                d[item] = kls(l.text)
        #print(d)
        return cls(**d)

    def toETree(self,et):
        loc=self.__class__.getSection(et)
        for fieldname,value in self.items():
            if isinstance(value,str):
                loc.find(fieldname).text=value
            elif isinstance(value,Number):
                loc.find(fieldname).text=str(value)
            elif isinstance(value,AbstractConfigSection):
                value.toETree(loc.find(fieldname))
            else:
                raise ValueError(str(fieldname) + " " +str(value))


class MicroBaseRates(AbstractConfigSection):
    __slots__ = ("decay","affinity","frequency",
    "prodEffectStrength","decayEffectStrength")
    __SlotClasses__=(ContentType.FLOAT,ContentType.FLOAT,ContentType.FLOAT,ContentType.FLOAT,
    ContentType.FLOAT)
    xmlSearchString="baseRates/micro"

    @classmethod
    def fromETree(cls,et):
        loc = cls.getSection(et)
        d={}
        for item,tp in zip(cls.__slots__,cls.__SlotClasses__):
            kls = clsCtorFor(tp)
            l=loc.find(item)
            #print("hello",kls,l.text)
            if issubclass(kls,AbstractConfigSection):
                d[item]=kls.fromETree(l)
            else:
                d[item] = kls(l.text)
        #print(d)
        return cls(**d)

    def toETree(self,et):
        loc=self.__class__.getSection(et)
        for fieldname,value in self.items():
            if isinstance(value,str):
                loc.find(fieldname).text=value
            elif isinstance(value,Number):
                loc.find(fieldname).text=str(value)
            elif isinstance(value,AbstractConfigSection):
                value.toETree(loc.find(fieldname))
            else:
                raise ValueError(str(fieldname) + " " +str(value))

class ConnectionParameters(AbstractConfigSection):
    __slots__ = ("miRNAs","TFs_coding","TFs_noncoding")
    __SlotClasses__ = (ContentType.DISCRETEDISTRIBUTION,ContentType.DISCRETEDISTRIBUTION,ContentType.DISCRETEDISTRIBUTION)
    xmlSearchString = "connectionParameters"
    
    @classmethod
    def getXMString(cls):
        xmlstrs=[item.replace("_","/") for item in cls.__slots__]
        return zip(cls.__slots__,xmlstrs)

    @classmethod
    def fromETree(cls,et):
        loc=cls.getSection(et)
        toRet = cls()
        for fieldname,xmlStr in cls.getXMString():
            #print(xmlStr)
            l=loc.find(xmlStr)
            #print(l.text)
            toRet[fieldname]=DiscreteDistributionXMLFactory(l)
        return toRet

    def toETree(self,et):
        loc=self.__class__.getSection(et)
        for fieldname,xmlStr in self.__class__.getXMString():
            l=loc.find(xmlStr)
            self[fieldname].toETree(l)

    @classmethod
    def generateConnectionsSweep(cls,key,**kwargs):
        if key not in cls.__slots__:
            raise AttributeError("key" +str(key) +" not in " +str(cls.__slots__))
        else:
            ind = cls.__slots__.index(key)
            sltype = clsFor()

class ConnectionType(Enum):
    TFCoding = "TFs/coding"
    TFNoncoding = "TFs/noncoding"
    miRNAs = "miRNAs"

    def xPathString(self):
        return "connectionParameters/"+self.value + "/discrete_distribution/distribution"

class KineticsVars(AbstractConfigSection):
    __slots__ = ("TXC__K_ONE_HALF","TXC__PROD_EXPONENT","TXL__K_ONE_HALF",
    "TXL__PROD_EXPONENT","MESS_DECAY__K_ONE_HALF","MESS_DECAY__EXPONENT",
    "MESSENGERDECAYFORCESMICRODECAY","MESSENGERDECAYFORCESMICROUNBINDING")
    __SlotClasses__ = (ContentType.FLOAT,ContentType.FLOAT,ContentType.FLOAT,
    ContentType.FLOAT,ContentType.FLOAT,ContentType.FLOAT,ContentType.BOOLEAN,ContentType.BOOLEAN)
    xmlSearchString="kinetics"

    @classmethod
    def getXMString(cls):
        xmlstrs=[item.replace("__","/") for item in cls.__slots__]
        return zip(cls.__slots__,xmlstrs)
    
    @classmethod
    def fromETree(cls,et):
        loc = cls.getSection(et)
        d={}
        for fieldname,xmlStr,tp in zip(*zip(*cls.getXMString()),cls.__SlotClasses__):
            kls = clsCtorFor(tp)
            l=loc.find(xmlStr)
            d[fieldname] = kls(l.attrib["value"])
        #print(d)
        return cls(**d)

    def toETree(self,et):
        loc=self.__class__.getSection(et)
        for fieldname,xmlStr in self.__class__.getXMString():
            value = self[fieldname]
            if isinstance(value,str):
                loc.find(xmlStr).attrib["value"]=value
            elif isinstance(value,Number):
                loc.find(xmlStr).attrib["value"]=str(value)
            elif isinstance(value,bool):
                loc.find(xmlStr).attrib["value"]=str(value).lower()
            elif isinstance(value,AbstractConfigSection):
                value.toETree(loc.find(xmlStr))
            else:
                raise ValueError(str(fieldname) + " " +str(value))


    def __setattr__(self,key,value):
        if value is None:
            super().__setattr__(key,value)
        elif key == "MESSENGERDECAYFORCESMICRODECAY":
            if not isinstance(value,bool):
                raise TypeError("value for MESSENGERDECAYFORCESMICRODECAY not bool")
            if value is True:
                super().__setattr__("MESSENGERDECAYFORCESMICRODECAY",True)
                super().__setattr__("MESSENGERDECAYFORCESMICROUNBINDING",False)
            else:
                super().__setattr__("MESSENGERDECAYFORCESMICRODECAY",False)
                super().__setattr__("MESSENGERDECAYFORCESMICROUNBINDING",True)
        elif key == "MESSENGERDECAYFORCESMICROUNBINDING":
            if not isinstance(value,bool):
                raise TypeError("value for MESSENGERDECAYFORCESMICROUNBINDING not bool")
            if value is False:
                super().__setattr__("MESSENGERDECAYFORCESMICRODECAY",True)
                super().__setattr__("MESSENGERDECAYFORCESMICROUNBINDING",False)
            else:
                super().__setattr__("MESSENGERDECAYFORCESMICRODECAY",False)
                super().__setattr__("MESSENGERDECAYFORCESMICROUNBINDING",True)
        else:
            super().__setattr__(key,value)



class ConfigFile:
    __slots__ = ("e","connections","kinetics","tfrates","messrates","microrates","ngs","modified")
    def __init__(self,fname):
        self.e = et.parse(fname)
        self.connections = ConnectionParameters.fromETree(self.e)
        self.kinetics = KineticsVars.fromETree(self.e)
        self.tfrates = TFBaseRates.fromETree(self.e)
        self.messrates = MessBaseRates.fromETree(self.e)
        self.microrates = MicroBaseRates.fromETree(self.e)
        self.ngs = self.e.find('.//randomSeeds/networkGenSeed')
        self.modified = False
    #def addSweep(self,key)

    #def scaleRatesToK(self):
        #km=self.kinetics.TXC__K_ONE_HALF
        #self.messrates.production *= 

    def toETree(self):
        self.connections.toETree(self.e)
        self.kinetics.toETree(self.e)
        self.tfrates.toETree(self.e)
        self.messrates.toETree(self.e)
        self.microrates.toETree(self.e)

    def renewNGS(self):
        netGenSeedVal = rnd.randint(0,100000)
        self.ngs.text = str(netGenSeedVal)
        self.modified=True

    def sweepTFEffectRates(self,mu=None,scale=None,sigma=None,displacement=None):
        argdict={}
        if mu is not None:
            argdict["mu"]=mu
        if scale is not None:
            argdict["scale"]=scale
        if sigma is not None:
            argdict["sigma"]=sigma
        if displacement is not None:
            argdict["displacement"]=displacement
        if len(argdict):
            return self.tfrates.effect.generateSweep(**argdict)

    def sweepNBConns(self,nSuccessRange=None,successPRange=None):
        argdict={}
        if nSuccessRange is not None:
            argdict["nSuccess"]= nSuccessRange
        if successPRange is not None:
            argdict["successP"]=successPRange
        if len(argdict):
            g1=self.connections.miRNAs.generateSweep(**argdict)
            g2=self.connections.TFs_coding.generateSweep(**argdict)
            g3=self.connections.TFs_noncoding.generateSweep(**argdict)
            return product(g1,g2,g3)
            


    def setMiRProdAndDecayEff(self,mPE,mDE):
        self.microrates.prodEffectStrength=mPE
        self.microrates.decayEffectStrength=mDE


    def setAndExec(self,nNetworks,argstup,mainstr=None):
        for item in argstup:
            #print(item[0],item[1])
            self.tfrates.effect.__setattr__(item[0],item[1])
        for i in range(nNetworks):
            self.execute(renewGenSeed=True,mainstr=mainstr)

    def write(self,fname):
        self.e.write(fname)
    

    def execute(self,renewGenSeed=True,mainstr=None):
        #if renewGenSeed:
        if mainstr is not None:
            mainnetgenstr = "./"+mainstr
        else:
            mainnetgenstr = "./mainnetgen"
        pargs = [mainnetgenstr,"-g"]
        self.renewNGS()
        if self.modified:
            self.toETree()
            self.write("config.xml")
            self.modified = False
        FNULL=open(os.devnull,'w')
        subprocess.call(pargs,stdout=FNULL,stderr=subprocess.STDOUT)

    #def 


def setConnection(e,g,vals):
    if not isinstance(e,et):
        raise TypeError("First argument must be XML Element tree")
    if not isinstance(g,ConnectionType):
        raise TypeError("Second argument must be ConnectionType")
    if not isinstance(vals,ConnectionDistribution):
        raise TypeError("Third argument must be ConnectionDistribution")
    distn=e.find(g.xPathString())
    for k,v in vals._asdict().items():
        if v is not None:
            distn.set(k,str(v))

def getConnection(e,g):
    if not isinstance(e,et):
        raise TypeError("First argument must be XML Element tree")
    if not isinstance(g,ConnectionType):
        raise TypeError("Third argument must be ConnectionType")
    distn=e.find(g.xPathString())
    return ConnectionDistribution(**distn.attrib)

def xmlnamesforKinVars(d=KineticsVars):
    return zip(xmlsearchnames,d)

def setTFEffectsLogNorm(e=et,vals=TFEffectsLogNorm):
    
    distn=e.find(tfeffectxmlsearchname)
    for k,v in vals._asdict().items():
        if v:
            distn.set(k,str(v))

def getTFEffectsLogNorm(e=et):
    distn = e.find(tfeffectxmlsearchname)
    d=dict(distn.attrib)
    for f in TFEffectsLogNorm._fields:
        if f not in d:
            d[f]=0
        elif not d[f].isalpha():
            d[f]=float(d[f])
    return TFEffectsLogNorm(**d)

def setKineticsVars(e=et,vals=KineticsVars): 
        y=xmlnamesforKinVars(vals)
        for k,v in y:
            if v is not None:
                e.find(k).set("value",str(v))

def getKineticsVars(e=et):
    j=[]
    for item in xmlsearchnames:
        j.append(e.find(item).get("value"))
    for n,item in enumerate(j):
        if item.lower() == "true":
            j[n]=True
        elif item.lower() == "false":
            j[n]=False
        else:
            j[n]=float(item)
    return KineticsVars(*j)



class cmdrun:
    mainfmt='mainsim.{0}.{1}.{2}.{3}.{4}.{5}.{6}'
    networkgenfmt='mainnetgen.{0}.{1}.{2}.{3}.{4}.{5}.{6}'
    builddirfmt='build.{0}.{1}.{2}.{3}.{4}.{5}.{6}'

    @staticmethod
    def formatMainNamesExt(kindict):
        #print(kindict)
        s=""
        if kindict["MESSENGERDECAYFORCESMICRODECAY"].lower()=="true":
            s="DEC"
        else:
            s="UNB"
        return cmdrun.mainfmt.format(kindict["TXC_K_ONE_HALF"],kindict["TXC_PROD_EXPONENT"],kindict["TXL_K_ONE_HALF"],kindict["TXL_PROD_EXPONENT"],
        kindict["MESS_DECAY_K_ONE_HALF"],kindict["MESS_DECAY_EXPONENT"],s)
        


        

    def __init__(self):
        self.cmdargs = {}

    def clearArgs(self):
        self.cmdargs = {}

    # def rebuild(self,withIPO=False):
    #     args=["make","clean"]
    #     subprocess.call(args)
    #     args=["make","opt"]
    #     if(withIPO):
    #         args.append('CMDLINEARGS="-ipo"')
    #     subprocess.call(args)

    def runSim(self,netgen=None):
        qargs=zip(self.cmdargs.keys(),self.cmdargs.values())
        if netgen:
            pargs=['./'+self.networkgen]
        else:
            pargs=['./'+self.main]
        pargs.extend(list(chain(*qargs)))
        FNULL=open(os.devnull,'w')
        def timing():
            #print(pargs)
            popen = subprocess.call(pargs,stdout=FNULL,stderr=subprocess.STDOUT)
            #popen.wait()
            return
        t=timeit.Timer(timing)
        time_taken=t.timeit(1)
        return time_taken

    def addArgs(self,flag,arg=None):
        self.cmdargs[flag]=str(arg)

    def setJSONFile(self,fname):
        if "genes.json" in fname:
            dname = fname[0:fname.find("genes.json")]
        else:
            dname=fname
        fn = os.path.join(dname,"config.xml")
        g=Configurator.kineticvars(fn)
        self.setkinetics(g)
        self.addArgs("-f",fname)

    def setNActiveMirs(self,n):
        self.addArgs("-n",n)

    def setNSimulations(self,n):
        self.addArgs("-s",n)

    def setGenerateNetworks(self):
        g=Configurator.kineticvars("config.xml")
        self.setkinetics(g)
        #print(self.networkgen)
        self.addArgs("-g")

    def formatMainNames(self,kindict):
        s=""
        if kindict["MESSENGERDECAYFORCESMICRODECAY"]=="true":
            s="DEC"
        else:
            s="UNB"
        self.networkgen=self.networkgenfmt.format(kindict["TXC_K_ONE_HALF"],kindict["TXC_PROD_EXPONENT"],kindict["TXL_K_ONE_HALF"],kindict["TXL_PROD_EXPONENT"],
        kindict["MESS_DECAY_K_ONE_HALF"],kindict["MESS_DECAY_EXPONENT"],s)
        self.main=self.mainfmt.format(kindict["TXC_K_ONE_HALF"],kindict["TXC_PROD_EXPONENT"],kindict["TXL_K_ONE_HALF"],kindict["TXL_PROD_EXPONENT"],
        kindict["MESS_DECAY_K_ONE_HALF"],kindict["MESS_DECAY_EXPONENT"],s)
        self.builddir=self.builddirfmt.format(kindict["TXC_K_ONE_HALF"],kindict["TXC_PROD_EXPONENT"],kindict["TXL_K_ONE_HALF"],kindict["TXL_PROD_EXPONENT"],
        kindict["MESS_DECAY_K_ONE_HALF"],kindict["MESS_DECAY_EXPONENT"],s)

    def setkinetics(self,kindict):
        self.formatMainNames(kindict)
        if not os.path.exists(self.networkgen):
            l=[]
            for k,v in kindict.items():
                if v=="true":
                    l.append("-D{0}".format(k))
                elif v=="false":
                    continue
                else:
                    l.append("-D{0}={1}".format(k,int(v)))
            defs=" ".join(l)
            
            #mkclean = subprocess.call(["make","clean"])
            parg='make CMDLINEARGS=\"' + defs + '\" TARGETSIM={0} TARGETNETGEN={1} BUILDDIR_ROOT={2}'.format(self.main,self.networkgen,
                self.builddir)
            #print(parg)
            popen = subprocess.call(parg,shell=True)
            shutil.rmtree(self.builddir)

def isDoneAlmost(dname,nAc,nSim,nTimes):
    if "ngs" not in dname:
        raise FileNotFoundError("ngs not in " + dname)
    elif "miRs" in dname:
        raise FileNotFoundError("mirs is in " + dname)
    mirdname="miRs.{0}".format(nAc)
    if mirdname not in os.listdir(dname):
        return None
    else:
        mirdname = os.path.join(dname,mirdname)
        nTimesDone = [z for z in os.listdir(mirdname) if z.startswith("states")]
        if nAc == 0:
            state0d=os.path.join(mirdname,nTimesDone[0])
            nSimsToDo= len([x for x in os.listdir(state0d) if x.endswith(".prot.txt")])
            nSimsNow = nSim - nSimsToDo
            return (nSimsNow,nTimes)
        nTimesNow= nTimes - len(nTimesDone)
        return (nSim,nTimesNow)

def timeSubprocessCall(args):
    timeStarted=time.time()
    FNULL=open(os.devnull,"w")
    subprocess.call(args,stdout=FNULL,stderr=subprocess.STDOUT)
    timeDelta=time.time()-timeStarted
    print("Time taken = {0}\t{1}".format(timeDelta,str(args)))

def doMain(args):
    fname=args[0] 
    nAc=args[1]
    nSim=args[2]
    nTimes=args[3]
    #print(args)
    if "genes.json" in fname:
        dname = fname[0:fname.find("genes.json")]
    else:
        dname=fname
    nSnT = isDoneAlmost(dname,nAc,nSim,nTimes)
    if nSnT is not None:
        if nSnT[0] == 0:
            return
        else:
            nSim=nSnT[0]
            nTimes=nSnT[1]
    fn = os.path.join(dname,"config.xml")
    try:
        g=Configurator.kineticvars(fn)
        qargs = [("-f",fname),("-n",str(nAc)),("-s",str(nSim))]
        mainname = cmdrun.formatMainNamesExt(g)
        pargs = ["./"+mainname]
        #FNULL = open(os.devnull,"w")
        pargs.extend(list(chain(*qargs)))
        if nAc == 0:
            timeSubprocessCall(pargs)
            #subprocess.call(pargs,stdout=FNULL,stderr=subprocess.STDOUT)
        else:
            for i in range(nTimes):
                timeSubprocessCall(pargs)
                #subprocess.call(pargs,stdout=FNULL,stderr=subprocess.STDOUT)
    except Exception as e:
        print("Error w/ "+ "\t".join(map(str,pargs)))
        print(e,str(e))

def remake(kindict,printRna=True):
    mainfmt='./mainsim.{0}.{1}.{2}.{3}.{4}.{5}.{6}'
    networkgenfmt='./mainnetgen.{0}.{1}.{2}.{3}.{4}.{5}.{6}'
    builddirfmt='build.{0}.{1}.{2}.{3}.{4}.{5}.{6}'
    def formatMainNames_(kindict):
        s=""
        if kindict["MESSENGERDECAYFORCESMICRODECAY"]=="true":
            s="DEC"
        else:
            s="UNB"
        networkgenname=networkgenfmt.format(kindict["TXC_K_ONE_HALF"],kindict["TXC_PROD_EXPONENT"],kindict["TXL_K_ONE_HALF"],kindict["TXL_PROD_EXPONENT"],
        kindict["MESS_DECAY_K_ONE_HALF"],kindict["MESS_DECAY_EXPONENT"],s)
        mainname=mainfmt.format(kindict["TXC_K_ONE_HALF"],kindict["TXC_PROD_EXPONENT"],kindict["TXL_K_ONE_HALF"],kindict["TXL_PROD_EXPONENT"],
        kindict["MESS_DECAY_K_ONE_HALF"],kindict["MESS_DECAY_EXPONENT"],s)
        builddirname=builddirfmt.format(kindict["TXC_K_ONE_HALF"],kindict["TXC_PROD_EXPONENT"],kindict["TXL_K_ONE_HALF"],kindict["TXL_PROD_EXPONENT"],
        kindict["MESS_DECAY_K_ONE_HALF"],kindict["MESS_DECAY_EXPONENT"],s)
        if os.path.exists(networkgenname):
            return None
        else:
            return {"networkgen" : networkgenname,"main" : mainname,"builddir_root":builddirname}
    
    def setDefs_(kindict,printrna=True):
        l=[]
        for k,v in kindict.items():
            if v=="true":
                l.append("-D{0}".format(k))
            elif v=="false":
                continue
            else:
                l.append("-D{0}={1}".format(k,int(v)))
        if printrna:
            l.append("-DPRINT_RNA")
        defs=" ".join(l)
        return defs
    d=formatMainNames_(kindict)
    parg=''
    if d is not None:
        defs=setDefs_(kindict,printrna=printRna)
        parg = 'make CMDLINEARGS=\"' + defs + '\" TARGETSIM={0} TARGETNETGEN={1} BUILDDIR_ROOT={2}'.format(d["main"],d["networkgen"],d["builddir_root"])
        print(parg)
        FNULL=open(os.devnull,"w")
        popen = subprocess.call(parg,shell=True)#,stdout=FNULL)
        FNULL.close()
        shutil.rmtree(d["builddir_root"])

class Configurator():
    def __init__(self,pth=None):
        if pth is None:
            self.outpath="config.xml"
            self.tree=et.parse(self.outpath)
            self.ro = False
        else:
            self.ro = True
            self.outpath=os.path.join(pth,"config.xml")
            self.tree=et.parse(self.outpath)
        self.e=self.tree.getroot()
        self.netGenSeed = self.e.find('.//networkGenSeed')
        self.nMess = self.e.find('.//nMess')
        self.nMicro = self.e.find('.//nMicro')
        self.modded = False

    def modification_function(func):
        def wrap(self,*args,**kwargs):
            retVal = func(self,*args,**kwargs)
            self.modded=True
        return wrap

    @modification_function
    def renewNGS(self):
        netGenSeedVal = rnd.randint(0,100000)
        self.netGenSeed.text = str(netGenSeedVal)        



    def setMiRProdAndDecayEff(self,mPE,mDE):
        #print("setmirprodandecy {0} {1}".format(mPE,mDE))
        self.mirProdEffectVal = mPE
        self.mirDecayEffectVal = mDE
        self.mirDecayEffect.text = str(self.mirDecayEffectVal)
        self.mirProdEffect.text = str(self.mirProdEffectVal)
        self.write()


    def setMicroProd(self,mp):
        self.microprodval=mp
        self.microproduction.text = str(self.microprodval)
        self.write()

    def getNMicroInterval(self):
        ival = 24//4
        if 24 % ival == 0:
            tmax = 25
        else:
            tmax = 25
        return list(range(0,tmax,ival))

    def write(self):
        if self.ro:
            raise ValueError("cfg in readonly mode")
        else:
            self.tree.write(self.outpath)
    
    def setMessProductionBaseTarget(self):
        pass
    
    @staticmethod
    def getKineticsVars_(etroot): 
        d={} 
        k=etroot.find("kinetics") 
        for kintype in k: 
            pfx=kintype.tag 
            if len(kintype): 
                for child in kintype: 
                    mykey="_".join((pfx,child.tag)) 
                    d[mykey]=child.attrib["value"] 
            else: 
                myval=kintype.get("value",None) 
                if myval: 
                    d[kintype.tag]=myval                       
        return d 

    def getKineticsVars(self):
            return Configurator.getKineticsVars_(self.e)

    @staticmethod
    def kineticvars(fname):
        tree = et.parse(fname)
        rt = tree.getroot()
        return Configurator.getKineticsVars_(rt)





runCommand = cmdrun()
cfg = Configurator(None)

def getDirs(origset=None,pth='./output'):
    snew= set([os.path.join(pth,x) for x in os.listdir(pth) if (os.path.isdir(os.path.join(pth,x)) & (x.startswith('net.md')))])
    if origset:
        return snew.difference(origset)
    else:
        return snew



def genNetwork(nNetworks,txck=None,txcexp=None,txlk=None,txlexp=None,deck=None,deckexp=None,messforcesdecay=True,):
    runCommand.clearArgs()
    dirsold=getDirs()
    cfg.setKineticsVars(txck,txcexp,txlk,txlexp,deck,deckexp,messforcesdecay)
    cfg.write()
#     print("txck={} txcexp={} txlk={},txlexp={},deck={},deckexp={},messforcesdecay={}".format( txck,txcexp,txlk
# ,txlexp
# ,deck
# ,deckexp
# ,messforcesdecay))
    #print(Configurator.kineticvars("config.xml"))
 
    for i in range(nNetworks):
        cfg.renewNGS()
        cfg.write()
        #print(cfg.ro)
        runCommand.setGenerateNetworks()
        runCommand.runSim(netgen=True)

        if i==10:
            break
    dirsnew=list(getDirs(dirsold))
    if len(dirsnew):
        shutil.copy('config.xml',dirsnew[0])
    return dirsnew

FNNULL = open(os.devnull,"w")


def genNetworkWmain(nNetworks,netgenstr,txck=None,txcexp=None,txlk=None,txlexp=None,deck=None,deckexp=None,messforcesdecay=True,
mu=None,scale=None,sigma=None,displacement=None):
    #runCommand.clearArgs()
    dirsold=getDirs()
    cfg.setKineticsVars(txck,txcexp,txlk,txlexp,deck,deckexp,messforcesdecay)
    cfg.setTFEffects(mu,scale,sigma,displacement)
    cfg.write()
#     print("txck={} txcexp={} txlk={},txlexp={},deck={},deckexp={},messforcesdecay={}".format( txck,txcexp,txlk
# ,txlexp
# ,deck
# ,deckexp
# ,messforcesdecay))
    #print(Configurator.kineticvars("config.xml"))
    pargs = ['./'+netgenstr,'-g']
    for i in range(nNetworks):
        cfg.renewNGS()
        cfg.write()
        #print(cfg.ro)
        #runCommand.setGenerateNetworks()
        #runCommand.runSim(netgen=True)
        subprocess.call(pargs,stdout=FNNULL)

    dirsnew=list(getDirs(dirsold))
    if len(dirsnew):
        shutil.copy('config.xml',dirsnew[0])
    return dirsnew

def getNetworks(dirsnew):
    ngss = []
    print(dirsnew)
    #netsByP =[os.path.join('output',x) for x in os.listdir('output') if x.startswith('net.md')]
    for net in dirsnew:
        #print net
        #print os.listdir(net)
        ngss.extend([os.path.join(net,z) for z in os.listdir(net) if z.startswith('ngs')])
    #subNets = [[os.path.join(q,z) for z in os.listdir(q) if z.startswith('ngs.')] for q in netsByP]
    return ngss




import time
nProc=cpu_count()-4
rcs = [cmdrun()]*nProc



def runWorker(args):
    pth=args[0]
    nSims=args[1]
    nMirs=args[2]
    rc=args[3]
    #mylen=args[4]
    #print(args)
    #print("{0} of {1} nets".format(ctr.value,mylen))
    with ctr.get_lock():
        ctr.value +=1
    if ctr // 10 == 0:
        print(ctr)
    #sl = rnd.random()
    #print 'sleepstart - {0}'.format(sl)
    #time.sleep(sl*10)
    rc.clearArgs()
    rc.setJSONFile(pth)
    rc.setNActiveMirs(nMirs)
    rc.setNSimulations(nSims)
    #print 'sleepend - {0}'.format(sl)
    tt=rc.runSim()
    #("runtime for {0} miRs = {1}s".format(nMirs,tt))




from itertools import product,repeat


def getMicroInt():
    return cfg.getNMicroInterval()

def runSimsForNetworks(nNetworks,nSims,nMess):
    #print "runSimsStart"
    #runCommand.rebuild(False)
    dnew=genNetwork(nNetworks)
    #print "genNetwork"
    netpaths=getNetworks(dnew)
    #runCommand.rebuild(True)
    ival=cfg.getNMicroInterval()
    arglist=product(netpaths,(nSims),ival)
    for i in range(len(netpaths)):
        netpath = netpaths[i]
        #print ival
        for nMicro in ival:
            #print netpath,nMicro
            runSimsForNMicro(netpath,nSims,nMicro)
    print("{0} of {1} networks complete".format(i,len(netpaths)))


def make_kindict(txck,txcexp,txlk,txlexp,deck,deckexp,forcesdecay):
    keys = ["TXC_K_ONE_HALF","TXC_PROD_EXPONENT","TXL_K_ONE_HALF","TXL_PROD_EXPONENT",
        "MESS_DECAY_K_ONE_HALF","MESS_DECAY_EXPONENT","MESSENGERDECAYFORCESMICRODECAY","MESSENGERDECAYFORCESMICROUNBINDING"]
    vals=[txck,txcexp,txlk,txlexp,deck,deckexp]
    if forcesdecay:
        vals.append("true")
        vals.append("false")
    else:
        vals.append("false")
        vals.append("true")
    return dict(zip(keys,vals))

def generateNetworkWorker(args):
    kd = make_kindict(*args[:-1])
    #print(kd)
    #set
    #print('{0} of {1} mains generated'.format(ctr.value,mylen))
    #with ctr.get_lock():
    #    ctr.value +=1
    remake(kd)

def parseMain(confile,mainstr,nNets,mPRange,mERange,murange,scalerange,sigmarange,disprange):
    j=mainstr.split('.')
    junb = j[-1]
    jivals = j[1:-1]
    jivals = list(map(int,jivals))
    if junb == "UNB":
        fd  = False
    elif junb == "DEC":
        fd = True
    else:
        raise ValueError("invalid tail")
    jivals.append(fd)
    q=list(product(mPRange,mERange))
    for item in q:
        confile.setMiRProdAndDecayEff(item[0],item[1])
        s=confile.sweepTFEffectRates(mu=murange,scale=scalerange,sigma=sigmarange,displacement=disprange)
        for elem in s:
            confile.setAndExec(2,elem,mainstr=mainstr)

def testKinetics(nNetworks,txckhlist=(None,),txcprodelist=(None,),txlkhlist=(None,),txlexplist=(None,),deckhlist=(None,),deckexplist=(None,),
forcesdecay=(True,False)):
    ptot=list(product(txckhlist,txcprodelist,txlkhlist,txlexplist,deckhlist,deckexplist,forcesdecay))
    myNets=set()
    #print(ptot)
    for pr in ptot:
        #print(pr)
        d=genNetwork(nNetworks,txck=pr[0],txcexp=pr[1],txlk=pr[2],txlexp=pr[3],deck=pr[4],deckexp=pr[5],messforcesdecay=pr[6])
        d=getNetworks(d)
        #print(myNets)
        myNets.update(d)
    #print(myNets)
    #runAllNetworksAtOnce(list(myNets),nSims)

def foreachnetgen():
    confile = ConfigFile("config.xml")
    mng = [x for x in os.listdir('.') if x.startswith("mainnetgen.")]
    myPRan=np.log(np.logspace(0.01,0.1,2))
    myDRan= np.linspace(1.1,2.5,2)
    #myPRan = (myPRan[1],)
    #myDRan = (myDRan[1],)
    mymu=np.linspace(0.5,3,2)
    myscale=np.linspace(0.25,0.75,2)
    mysigma=np.linspace(0.75,1.25,2)
    mydisplacement=np.linspace(0,0.4,2)

    ct=0
    nng = len(mng)
    for item in mng:
        parseMain(confile,item,2,myPRan,myDRan,mymu,myscale,mysigma,mydisplacement)
        ct+=1
        print("{0} of {1}".format(ct,nng))

if __name__ == "__main__":
    cfg.renewNGS()

    # for r,d,f in os.walk('output'):
    #     if "genes.json" in f:
    #         nets.append(r)
    # nets=list(set(nets))
    # print(nets)
    # runAllNetworksAtOnce(nets,2)
    #tfiqrange=[100,500,1000,5000,10000,100000,500000]
    # prodrange = [0.25, 0.75, 1.25, 1.75, 2.25]
    # myNets=set()
    # numNetsPerParam=5
    # for prodv in prodrange:
    #     print(prodv)
    #     cfg.setMicroProd(prodv)
    #     d=genNetwork(numNetsPerParam)
    #     d=getNetworks(d)
    #     myNets.update(d)
    # print(myNets)
    # runAllNetworksAtOnce(list(myNets),1)
