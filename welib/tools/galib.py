""" 
 Genetic algorithm tools
 Uses the same conventions as DEAP:
    fitness values are stored in 
        p.fitness.values
    p is a list
"""

import os
import numpy as np
import random
from copy import deepcopy
from collections import Sequence
from itertools import repeat
import hashlib
import math
import glob
import re
import pandas as pd

class Fitn():
    pass

class Indiv(list):
    def __init__(self, *args):
        list.__init__(self, *args)
        self.fitness=Fitn()
        self.fitness.values=[]
# --------------------------------------------------------------------------------}
# --- Gene MAP
# --------------------------------------------------------------------------------{
class GeneMap():
    def __init__(self,name,nBases,protein_ranges,kind=None,protein_neutr=None,meta=None,resolution=1000):
        """
        A gene is between 0 and 1 
        A protein is defined by the ranges
        """
        self.nBases = nBases
        self.kind   = kind
        self.protein_ranges = protein_ranges
        self.protein_neutr  = protein_neutr
        if protein_neutr is None:
            self.protein_neutr=[(m+M)/2 for m,M in protein_ranges]
        self.meta  = meta
        self.resolution  = resolution

        def pretty_name(n):
            if n.find('|')>0:
                s=n.split('|')
                return s[-1]
            elif n.find('\\')>0:
                n=n.replace('.dat','')
                s=n.split('\\')
                return s[-1]
            elif n.find('/')>0:
                n=n.replace('.dat','')
                s=n.split('/')
                return s[-1]
            else:
                return n
        self.name          = name
        self.pretty_name   = pretty_name(name)

    def __repr__(self):
        s=''.join(['x']*self.nBases)+': '+self.name+'\n'
        return s

    def decode(self,gene,iBase=None):
        if iBase is None:
            prot =[]
            for g,pr in zip(gene,self.protein_ranges):
                p=pr[0]+ g*(pr[1]-pr[0])
                prot.append(p)
                if g<0 or g>1:
                    print('g:',g, 'pr:',pr, ' -> p:',p)
                    raise Exception('The gene cannot be decoded properly')
        else:
            g=gene
            pr=self.protein_ranges[iBase]
            prot=pr[0]+ g*(pr[1]-pr[0])
            if g<0 or g>1:
                print('g:',g, 'pr:',pr, ' -> p:',prot)
                raise Exception('The base cannot be decoded properly')
        return prot

    def encode(self,protein):
        gene=[]
        for p,pr in zip(protein,self.protein_ranges):
            g=(p-pr[0])/(pr[1]-pr[0])
            gene.append(g)
            if g>1 or g<0:
                print('p:',p, 'pr:',pr, ' -> g:',g)
                raise Exception('The protein cannot be encoded properly')
        return gene

    def neutralProtein(self):
        return self.protein_neutr

    def show_full_raw(self,gene):
        s='['+' '.join([str(b) for b in gene])+']: '+self.name
        return s

    def show_full(self,gene):
        def pretty(b,pr):
            if self.resolution is None:
                return str(b)
            delta=pr[1]-pr[0]
            if delta<=0:
                return str(b)
            nDec=int(np.log10(delta/self.resolution))
            nInt=int(np.log10(pr[1]))
            if nInt<0:
                nInt=-1
            if nDec<0:
                fmt='{:'+str(nInt-nDec+3)+'.'+str(-nDec+1)+'f}'
                #print(fmt)
                return fmt.format(b)
            elif nInt>0:
                fmt='{:'+str(nInt+1)+'.0f}'
                #print(fmt)
                return fmt.format(b)
            else:
                return str(b)
        
        if self.nBases>1:
            s=self.pretty_name+': ['+' '.join([pretty(b,rg) for b,rg in zip(self.decode(gene),self.protein_ranges)])+']'
        else:
            s=self.pretty_name+': '+pretty(self.decode(gene)[0],self.protein_ranges[0])
        return s

    def geneBounds(self):
        return [(0,1)]*self.nBases

    def proteinBounds(self):
        return [(m,M) for m,M in self.protein_ranges]



class ChromosomeMap(list):
    def add(self,d):
        self.append(d)

    def append(self,d):
        super(ChromosomeMap,self).append(d)
        # TODO check that 
        if not isinstance(d,GeneMap):
            raise Exception('Can only add `GenMap` types')

    @property
    def nBases(self):
        return sum([gene.nBases for gene in self])

    @property
    def nGenes(self):
        return len(self)

    def maxNBasesPerGene(self):
        return max([gene.nBases for gene in self])

    def neutralChromosome(self):
        v=[]
        for gm in self:
            v+=gm.encode(gm.protein_neutr)
        return v

    def neutralProtein(self):
        v=[]
        for gm in self:
            v+=gm.protein_neutr
        return v

    def decode(self,chromosome,iBase=None):
        if iBase is None:
            v=[]
            for gm,gene in zip(self,self.split(chromosome)):
                v+=gm.decode(gene) 
        else:
            if iBase>=self.nBases:
                raise Exception('iBase should be between 0 and nBases')
            i=0
            for ig,gm in enumerate(self):
                if (iBase>=i) and (iBase<i+gm.nBases):
                    break
                else:
                    i+=gm.nBases
            iBase=iBase-i
            #print('New iBase: {} for gene: {} {}'.format(iBase,ig,gm))
            v=gm.decode(chromosome,iBase)
        return v

    def encode(self,protein_chain):
        v=[]
        for gm,prot in zip(self,self.split(protein_chain)):
            v+=gm.encode(prot) 
        return v

    def chromosomeBounds(self):
        v=[]
        for gm in self:
            v+=gm.geneBounds() 
        return v

    def proteinChainBounds(self):
        v=[]
        for gm in self:
            v+=gm.proteinBounds() 
        return v


    def __repr__(self):
        s=''
        fmt='{:'+str(self.maxNBasesPerGene())+'s}'
        for gm in self:
            s+=fmt.format(''.join(['x']*gm.nBases))+': '+gm.name+'\n'
        return s

    def show_full(self,chromosome,sep='\n'):
        s=''
        for gm,gene in zip(self,self.split(chromosome)):
            s+=gm.show_full(gene)+sep
        return s

    def split(self,chromosome):
        genes=[]
        n=0
        for gm in self:
            genes.append(chromosome[n:n+gm.nBases])
            n=n+gm.nBases
        return genes

# --------------------------------------------------------------------------------}
# --- ID  
# --------------------------------------------------------------------------------{
def nparray_hash(x,length=16):
   return hashlib.md5((x.tobytes())).hexdigest()[:length]

def chromID(p):
    return nparray_hash(np.array(p),length=32)

# --------------------------------------------------------------------------------}
# --- Parametric 
# --------------------------------------------------------------------------------{
def parameticGA(fitnessEvalFun,ch_map,nPerBase,nFitness,resolution=None):
    """ 
        Perform a parametric study using the same formalism of the Genetic algorithm
        Each base is varies between 0 and 1 as defined by `nPerBase` (a list of values for each base or a single value)
        The function `fitnessEvalFun` is evaluated on the population

        `resolution` should be a power of 10, like 10, 100, 1000
    
    """
    nBases=ch_map.nBases
    if isinstance(nPerBase,list):
        if len(nPerBase)!=nBases:
            raise Exception('If nPerBase is a list it must be the same length as the number of bases')
    else:
        nPerBase= [nPerBase]*nBases

    nTot       = np.prod(nPerBase)
    nValuesCum = np.insert(np.cumprod(nPerBase),0,1)[:-1];
    vBaseValues=[np.linspace(0,1,n) for n in nPerBase]
    print('Parametric values (no rounding:)')
    for v in vBaseValues:
        print(v)
    vProtValues=[np.array([ch_map.decode(g,iBase=j) for g in v]) for j,v in enumerate(vBaseValues)]
    print('Prot values (no rounding:)')
    for v in vProtValues:
        print(v)
    if resolution:
        vBaseValues=[np.round(resolution*np.linspace(0,1,n))/resolution for n in nPerBase]
        print('Parametric values (with rounding:)')
        for v in vBaseValues:
            print(v)
        # we scale 
        print('Prot values (with rounding:)')
        vProtValues=[np.array([ch_map.decode(g,iBase=j) for g in v]) for j,v in enumerate(vBaseValues)]
        for v in vProtValues:
            print(v)

    fits_arr  = np.zeros( tuple(nPerBase+[nFitness] ) )
    fits_norm = np.zeros( tuple(nPerBase) )
    ValFit    = np.zeros( tuple(nPerBase) )
    
    print('Creating population of {} individuals...'.format(nTot))
    pop=[]
    for i in range(nTot):
        Indexes=(np.mod(np.floor(i/nValuesCum),nPerBase)).astype(int);
        chromosome=Indiv([vBaseValues[j][Indexes[j]] for j in range(nBases) ])
        #print(i,Indexes,chromosome)
        pop.append(chromosome)

    print('Evaluating population...')
    for i,p in enumerate(pop):
        Indexes=tuple((np.mod(np.floor(i/nValuesCum),nPerBase)).astype(int));
        fits = fitnessEvalFun(p,stat='{:4.1f}% - '.format(100.0*i/nTot))
        fits_norm[Indexes] = np.linalg.norm(fits)
        fits_arr [Indexes] = fits
    return fits_norm,fits_arr,pop,vBaseValues,vProtValues

# --------------------------------------------------------------------------------}
# --- Population/individual manipulations
# --------------------------------------------------------------------------------{
def addIfAbsent(pop,ind,sPop='',sInd=''):
    if ind not in pop:
        pop.append(clone(ind))
        if len(sPop)>0:
            print('Adding {:8s} to {}'.format(chromID(ind),sPop))

def clone(x):
    return deepcopy(x)


def populationTrimAccuracy(pop,nDecimals=None):
    for i in range(len(pop)):
        for j in range(len(pop[0])):
            pop[i][j]=np.around(pop[i][j],decimals=nDecimals)
    return pop


def splitstrn(s,n):
    return [s[i:i+n] for i in range(0, len(s), n)]


def populationStats(pop,best=None,stats=None):
    nVar = len(pop[0].fitness.values)
    new_stats=[]
    for i in range(nVar):
        fits = [p.fitness.values[i] for p in pop]
        d=dict()
        d['Mean'] = np.mean(fits)
        d['Min']  = np.min(fits)
        d['Max']  = np.max(fits)
        d['Std']  = np.std(fits)
        if best is not None:
            d['Best'] = best.fitness.values[i]
        else:
            d['Best'] = np.nan
        new_stats.append(pd.DataFrame(d,index=[0]))
    if stats is not None:
        stats=[ s.append(ns, ignore_index=True) for s,ns in zip(stats,new_stats)]

    return new_stats,stats

def populationPrint(pop,nBasePerGene=None,label=''):
    print('------------------ {} POPULATION ----------------------'.format(label))
    if nBasePerGene is None:
        nBasePerGene=len(pop[0])
    for p,i in zip(pop,range(len(pop))):
        nCharPerBase = 5 # NOTE related to format below..
        sBases=' '.join(['{:.2f}'.format(x) for x in p])
        splits=splitstrn(sBases,nBasePerGene*nCharPerBase)
        sGenes='| '.join(splits)
        sFits=' '.join(['{:.3f}'.format(x) for x in p.fitness.values])
        print('#{:2d} | {} | {}'.format(i,sFits,sGenes))

def populationSave(pop,directory='',basename='GA_DB_',newfile=False,fileformat='csv',fitsnames=None,basesnames=None,fitsformats=None,basesformats=None):
    # Detecting existing files
    files=glob.glob(os.path.join(directory,basename+'[0-9]*'+'.'+fileformat))
    if newfile:
        if len(files)==0:
            i=0
        else:
            i=int(re.search(r'\d+',files[-1]).group())+1
        filename=os.path.join(directory,'{}{:04d}.{}'.format(basename,i,fileformat))
        readflag='w'
    else:
        filename=files[-1]
        readflag='a'

    directory=os.path.dirname(filename)
    if len(directory)>0:
        if not os.path.exists(directory):
            os.mkdir(directory)

    # Creating a matrix with the population data  
    nRow   = len(pop)
    nBases = len(pop[0])
    nFit   = len(pop[0].fitness.values)

    if fitsformats is None:
        fitsformats='{:.5f}'
    if not isinstance(fitsformats, list):
        fitsformats=[fitsformats]*nFit
    if basesformats is None:
        basesformats='{:.5f}'
    if not isinstance(basesformats, list):
        basesformats=[basesformats]*nBases

    # writing to disk
    with open(filename,readflag) as f:
        if fileformat.lower()=='csv':
            delim=' '
            if newfile:
                # Header
                if fitsnames is None:
                    fitsnames=['Fit{:d}'.format(i) for i in range(nFit)]
                if basesnames is None:
                    basesnames=['Base{:d}'.format(i) for i in range(nBases)]
                s = 'ID'+delim+delim.join(fitsnames)+delim+delim.join(basesnames)
                f.write(s)
                f.write('\n')
            for i,p in enumerate(pop):
                sFits  = delim.join([s.format(v) for s,v in zip(fitsformats,p.fitness.values)])
                sBases = delim.join([s.format(v) for s,v in zip(basesformats,p)])
                f.write(chromID(p)+delim+sFits+delim+sBases+'\n')
        else:
            raise Exception('Unknown fileformat {}'.format(fileformat))
    return filename

def populationLoad(filename=None, nFits=2):
    fileformat=os.path.splitext(filename)[1][1:].lower()

    if fileformat == 'csv':
        df=pd.read_csv(filename, sep=' ')
    else:
        raise Exception('Unknown file format: {}'.format(fileformat))

    # Unique values
    dfu= df.drop_duplicates().copy()
    dfu.reset_index(drop=True,inplace=True);
    pop_unique  = dfu.values[:,(nFits+1):].astype(float)
    fits_unique = dfu.values[:,1:(nFits+1)].astype(float)

    class fit():
        pass
    class ind(list):
        def __init__(self,chromosome,fits):
            super(ind,self).__init__(chromosome)
            self.fitness=fit()
            self.fitness.values=fits
    pop=[ind(c,f) for c,f in zip(pop_unique,fits_unique)]
    #return df,pop_unique,fits_unique
    return df,pop


def timelinePlot(df,fig,noFit=True):
    if fig is None:
        fig=plt.figure(figsize=(5, 12));
    if len(fig.axes)==0:
        ax=fig.add_subplot(111);
    ax=fig.axes[0]
    ax.clear()
    # shifting values for plottting
    numeric_cols = [col for col in df if df[col].dtype.kind != 'O']
    if noFit:
        numeric_cols = [col for col in numeric_cols if col.find('Fit')!=0]
    for i,c in enumerate(numeric_cols):
        df[c]= df[c] + len(numeric_cols)-i-1
    df[numeric_cols].plot(ax=ax)
    for i in range(len(numeric_cols)):
        ax.plot([0,len(df)],[i,i],'k--');
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5));
    ax.set_title('Evolution of chromosome bases')
    ax.set_xlabel('Cumulative number of chromosomes')
    ax.set_ylabel('Bases number')



# --------------------------------------------------------------------------------}
# --- Selection
# --------------------------------------------------------------------------------{
def selBestByNorm(pop,kind=None,weights=None):
    if kind is None:
        # we determine the kind based on the sign of weights
        # minimization problem, negative weights => min of norm
        # maximization problem, positive weights => max of norm
        nVar=len(pop[0].fitness.values)
        if weights is None:
            weights=[1.0]*nVar
        signs=np.unique(np.sign(weights))
        if len(signs)>1:
            raise Exception('Select by norm can only work for all minimize or all maximize')
        sign=signs[0]
        if sign>0:
            kind='max'
        else: 
            kind='min'
    norms=[np.linalg.norm(np.array(p.fitness.values)) for p in pop]
    if kind=='max':
        i=np.argmax(norms)
    elif kind=='min':
        i=np.argmin(norms)
    else:
        raise Exception('Wrong kind')

    return clone(pop[i])


def selIndependentBest(pop,kind='max',weights=None):
    if weights is None:
        weights=[1.0]*len(pop[0].fitness.values)

    best_inds=[]
    for i in range(len(pop[0].fitness.values)):
        fits=[weights[i]*p.fitness.values[i] for p in pop]
        if kind=='max':
            i=np.argmax(fits)
        elif kind=='min':
            i=np.argmin(fits)
        else:
            raise Exception('Wrong kind')
        best_inds.append(clone(pop[i]))
    return best_inds
# --------------------------------------------------------------------------------}
# --- Mating 
# --------------------------------------------------------------------------------{
def mate(population,fMate,nChildren=None):
    # Cross over methods :
    #cxOnePoint() #cxTwoPoint() #cxUniform() #cxPartialyMatched() #cxUniformPartialyMatched() 
    #cxOrdered() #cxBlend() #cxESBlend() #cxESTwoPoint() #cxSimulatedBinary()
    #cxSimulatedBinaryBounded() #cxMessyOnePoint()
    nParents=len(population)
    if nChildren is None:
        nChildren=nParents

    children=[]
    for i in range(nChildren):
        dad = clone(population[random.randint(0,nParents-1)])
        mom = clone(population[random.randint(0,nParents-1)])
        child1,child2=fMate(dad, mom);
        if random.random()>0.5:
            child=child1
        else:
            child=child2
        if hasattr(child,'data'):
            del child.data
        del child.fitness.values
        children.append(child)
    return children


def mateInPlace(population,fMate,CXPB):
    # --- CROSSOVER of populations, with probability CXPB
    for child1, child2 in zip(population[::2], population[1::2]):
        if random.random() < CXPB:
            fMate(child1, child2);
            del child1.fitness.values
            del child2.fitness.values
            if hasattr(child1,'data'):
                del child1.data
            if hasattr(child2,'data'):
                del child2.data
    nCX = sum(1 for ind in population if not ind.fitness.valid)

# --------------------------------------------------------------------------------}
# --- Mutation 
# --------------------------------------------------------------------------------{
def mutateInPlace(pop,fMutate,MUTPB):
    nMute=0;
    for mutant in pop:
        if random.random() < MUTPB:
            fMutate(mutant);
            del mutant.fitness.values
            if hasattr(mutant,'data'):
                del mutant.data
            nMute+=1
    return nMute



def mutGaussianBounded(individual, mu, sigma, low, up, indpb):
    """This function applies a gaussian mutation of mean *mu* and standard
    deviation *sigma* on the input individual. This mutation expects a
    :term:`sequence` individual composed of real valued attributes.
    The *indpb* argument is the probability of each attribute to be mutated.
    :param individual: Individual to be mutated.
    :param mu: Mean or :term:`python:sequence` of means for the
               gaussian addition mutation.
    :param sigma: Standard deviation or :term:`python:sequence` of
                  standard deviations for the gaussian addition mutation.
    :param indpb: Independent probability for each attribute to be mutated.
    :returns: A tuple of one individual.
    This function uses the :func:`~random.random` and :func:`~random.gauss`
    functions from the python base :mod:`random` module.
    """
    size = len(individual)
    if not isinstance(mu, Sequence):
        mu = repeat(mu, size)
    elif len(mu) < size:
        raise IndexError("mu must be at least the size of individual: %d < %d" % (len(mu), size))
    if not isinstance(sigma, Sequence):
        sigma = repeat(sigma, size)
    elif len(sigma) < size:
        raise IndexError("sigma must be at least the size of individual: %d < %d" % (len(sigma), size))
    if not isinstance(low, Sequence):
        low = repeat(low, size)
    elif len(low) < size:
        raise IndexError("low must be at least the size of individual: %d < %d" % (len(low), size))
    if not isinstance(up, Sequence):
        up = repeat(up, size)
    elif len(up) < size:
        raise IndexError("up must be at least the size of individual: %d < %d" % (len(up), size))

    for i, m, s, xl, xu in zip(range(size), mu, sigma, low, up):
        if random.random() < indpb:
            dx=random.gauss(m, s)
            x=individual[i] + dx
            # if outside the limits, we bounce modulo the range
            if x>xu:
                x=xu-math.fmod(dx,xu-xl)
            elif x<xl:
                x=xl+math.fmod(abs(dx),xu-xl)
            #individual[i] = min(max(x,low),up)
            individual[i] = x

    return individual,

def mutUniformBounded(individual, low, up, indpb):
    """Mutate an individual by replacing attributes, with probability *indpb*,
    by a integer uniformly drawn between *low* and *up* inclusively.
    :param individual: :term:`Sequence <sequence>` individual to be mutated.
    :param low: The lower bound or a :term:`python:sequence` of
                of lower bounds of the range from wich to draw the new
                integer.
    :param up: The upper bound or a :term:`python:sequence` of
               of upper bounds of the range from wich to draw the new
               integer.
    :param indpb: Independent probability for each attribute to be mutated.
    :returns: A tuple of one individual.
    """
    size = len(individual)
    if not isinstance(low, Sequence):
        low = repeat(low, size)
    elif len(low) < size:
        raise IndexError("low must be at least the size of individual: %d < %d" % (len(low), size))
    if not isinstance(up, Sequence):
        up = repeat(up, size)
    elif len(up) < size:
        raise IndexError("up must be at least the size of individual: %d < %d" % (len(up), size))
    for i, xl, xu in zip(range(size), low, up):
        if random.random() < indpb:
            individual[i]= random.random()*(xu-xl)+xl


    return individual,


def chromosomeDistance(chrom1,chrom2):
    return np.linalg.norm(np.array(chrom1)-np.array(chrom2))

def populationDiversity(pop):
    return [sum([individualDistance(x,y) for x in pop])/(np.sqrt(2)*(len(pop)-1)) for y in pop]

def populationChromosomeDistances(pop):
    N = len(pop)
    distances = np.zeros((N,N))
    for i in range(N):
        for j in range(i+1, N):
            distances[i,j] = chromosomeDistance(pop[i],pop[j])
            distances[j,i] = distances[i,j]
    return distances


def selDiverse(pop, k):
    nInd= len(pop)
    nFit= len(pop[0].fitness.values)
    if nFit>1:
        raise Exception('Not intended for multi fitness')

    distances = populationChromosomeDistances(pop)
    distances = distances/np.max(distances)
    dist      = np.zeros(nInd)
    for i in range(nInd):
        dist[i] = np.sum(distances[i,:])

    dist = (dist-np.min(dist))/(np.max(dist)-np.min(dist))


    fits=np.array([p.fitness.values[0] for p in pop])
    fits = (fits-np.min(fits))/(np.max(fits)-np.min(fits))
    dist = dist*fits
    dist = (dist-np.min(dist))/(np.max(dist)-np.min(dist))

    fitA  = [(fits[i], i) for i in range(nInd)]
    distA = [(dist[i], i) for i in range(nInd)]
    fitA.sort(  reverse=True )
    distA.sort( reverse=True )

    #fitB=fits+fits*dist
    #fitBB = [(fitB[i], i) for i in range(nInd)]
    #fitBB.sort(  reverse=True )
    #BestI = [i for _,i in fitBB[:k] ]


    kBest=k-3
    if kBest<0:
        kBest=k
    BestFitI = [i for _,i in fitA[:kBest] ]
    BestDist = [i for _,i in distA if not i in BestFitI]
    BestDist = BestDist[:k-len(BestFitI)]
    nMissing = k-len(BestFitI)-len(BestDist)
    RandI = [random.randint(0,nInd) for i in range(nMissing)]
    BestI = BestFitI + BestDist + RandI
    if len(BestI)!=k:
        raise Exception('Screwed up')
    return [pop[i] for i in BestI]

def selSPEA2Diverse(individuals, k):
    """Apply SPEA-II selection operator on the *individuals*. Usually, the
    size of *individuals* will be larger than *n* because any individual
    present in *individuals* will appear in the returned list at most once.
    Having the size of *individuals* equals to *n* will have no effect other
    than sorting the population according to a strength Pareto scheme. The
    list returned contains references to the input *individuals*. For more
    details on the SPEA-II operator see [Zitzler2001]_.
    :param individuals: A list of individuals to select from.
    :param k: The number of individuals to select.
    :returns: A list of selected individuals.
    .. [Zitzler2001] Zitzler, Laumanns and Thiele, "SPEA 2: Improving the
       strength Pareto evolutionary algorithm", 2001.
    """
    N = len(individuals)
    nGenes= len(individuals[0])
    L = len(individuals[0].fitness.values)
    K = math.sqrt(N)
    strength_fits = [0] * N
    fits = [0] * N
    dominating_inds = [list() for i in range(N)]

    for i, ind_i in enumerate(individuals):
        for j, ind_j in enumerate(individuals[i+1:], i+1):
            if ind_i.fitness.dominates(ind_j.fitness):
                strength_fits[i] += 1
                dominating_inds[j].append(i)
            elif ind_j.fitness.dominates(ind_i.fitness):
                strength_fits[j] += 1
                dominating_inds[i].append(j)

    for i in range(N):
        for j in dominating_inds[i]:
            fits[i] += strength_fits[j]

    # Choose all non-dominated individuals
    chosen_indices = [i for i in range(N) if fits[i] < 1]

    if len(chosen_indices) < k:     # The archive is too small
        print('>>>>>> TOO SMALL', len(chosen_indices),k)
        distances = populationChromosomeDistances(individuals)
        distances=distances/np.max(distances)
        #[print('Chosen',chosen_indices)
        #[print('Ind',i)
        for i in range(N):
            print(distances[i,:])
            kth_dist = _randomizedSelect(distances[i,:], 0, N - 1, K)
            density = 1.0 / (kth_dist + 2.0)
            fits[i] += density

        next_indices = [(fits[i], i) for i in range(N) if not i in chosen_indices]
        next_indices.sort()
        #print next_indices
        chosen_indices += [i for _, i in next_indices[:k - len(chosen_indices)]]

    elif len(chosen_indices) > k:   # The archive is too large
        print('>>>>>> TOO BIG')
        N = len(chosen_indices)
        distances = [[0.0] * N for i in range(N)]
        sorted_indices = [[0] * N for i in range(N)]
        for i in range(N):
            for j in range(i + 1, N):
                dist = 0.0
                for l in range(L):
                    val = individuals[chosen_indices[i]].fitness.values[l] - \
                          individuals[chosen_indices[j]].fitness.values[l]
                    dist += val * val
                distances[i][j] = dist
                distances[j][i] = dist
            distances[i][i] = -1

        # Insert sort is faster than quick sort for short arrays
        for i in range(N):
            for j in range(1, N):
                l = j
                while l > 0 and distances[i][j] < distances[i][sorted_indices[i][l - 1]]:
                    sorted_indices[i][l] = sorted_indices[i][l - 1]
                    l -= 1
                sorted_indices[i][l] = j

        size = N
        to_remove = []
        while size > k:
            # Search for minimal distance
            min_pos = 0
            for i in range(1, N):
                for j in range(1, size):
                    dist_i_sorted_j = distances[i][sorted_indices[i][j]]
                    dist_min_sorted_j = distances[min_pos][sorted_indices[min_pos][j]]

                    if dist_i_sorted_j < dist_min_sorted_j:
                        min_pos = i
                        break
                    elif dist_i_sorted_j > dist_min_sorted_j:
                        break

            # Remove minimal distance from sorted_indices
            for i in range(N):
                distances[i][min_pos] = float("inf")
                distances[min_pos][i] = float("inf")

                for j in range(1, size - 1):
                    if sorted_indices[i][j] == min_pos:
                        sorted_indices[i][j] = sorted_indices[i][j + 1]
                        sorted_indices[i][j + 1] = min_pos

            # Remove corresponding individual from chosen_indices
            to_remove.append(min_pos)
            size -= 1

        for index in reversed(sorted(to_remove)):
            del chosen_indices[index]

    print(chosen_indices) 
    Sel=[individuals[i] for i in chosen_indices]
    print(len(chosen_indices),k)
    SelU=[]
    for i in chosen_indices:
        if individuals[i] not in SelU:
            SelU.append(individuals[i])

    print('Selected')
    print(len(Sel),k)
    #jjprint(Sel)
    print('Unique ones')
    print(len(SelU),k)
    #print(SelU)
    if len(SelU)<k:
        print('>>>>>> NEED FOR MORE')
        raise Exception()

    return  Sel



def _randomizedPartition(array, begin, end):
    i = random.randint(begin, end)
    array[begin], array[i] = array[i], array[begin]
    return _partition(array, begin, end)

def _partition(array, begin, end):
    x = array[begin]
    i = begin - 1
    j = end + 1
    while True:
        j -= 1
        while array[j] > x:
            j -= 1
        i += 1
        while array[i] < x:
            i += 1
        if i < j:
            array[i], array[j] = array[j], array[i]
        else:
            return j

def _randomizedSelect(array, begin, end, i):
    """Allows to select the ith smallest element from array without sorting it.
    Runtime is expected to be O(n).
    """
    if begin == end:
        return array[begin]
    q = _randomizedPartition(array, begin, end)
    k = q - begin + 1
    if i < k:
        return _randomizedSelect(array, begin, q, i)
    else:
        return _randomizedSelect(array, q + 1, end, i - k)
