import random
import os
import math
from itertools import combinations
import statistics
import rPAC.rpac_route as route 
import rPAC.rpac_score as score
import rPAC.rpac_pathway_graph as pathway_graph
import rPAC.rpac_util as util
from scipy.stats import pearsonr
import random


def singlecellmarkerlabel(namedsamples,constrain):

    '''
    Get qualified sample from namedsamples.
    Input : all named samples
            Constrain[genename][larger/smaller]=value
    '''
    selectedsampels={}
    for samplename,genedict in namedsamples.items():
        isadd=True
        for genename, contraindetail in constrain.items():
            for condition,value in contraindetail.items():
                if condition=='larger' and genedict[genename]<value:
                    isadd=False
                if condition=='smaller' and genedict[genename]>value:
                    isadd=False
        if isadd:
            selectedsampels[samplename]=genedict
    return selectedsampels

def loadfilewithname(inputfilepath,sep='\t' ,transpose=False):
    
    '''
    Load file from txt file

    First row should be sample name
    First column should be gene name
    other wise transpose should be True
    input 
        inputfilepath: Path to file
        transpose: Does this file need transpose. Default is transpose
    return a 2-D list. [name][genename] = float value 
    '''
    genelist = []
    allnamesamples={}
    with open(inputfilepath,'r') as allsamplefile:
        istitle=True
        if transpose:
            genelist = [] 
            for line in allsamplefile.readlines():
                linedata = line.strip().split(sep)
                if istitle:
                    del(linedata[0])
                    for genename in linedata:
                        genelist.append(genename.upper())                        
                    istitle=False
                else:
                    sampledetails ={} 
                    samplename =linedata[0] 
                    # print(samplename)
                    del(linedata[0])
                    for i in range(len(genelist)):                        
                        sampledetails[genelist[i]] = float(linedata[i])
                    allnamesamples[samplename]=sampledetails                    
        else:
            samplelist=[]

            for line in allsamplefile.readlines():
                linedata = line.strip().split(sep)
                if istitle:
                    del(linedata[0])
                    for sample in linedata:
                        allnamesamples[sample]={}
                        samplelist.append(sample)
                    istitle=False
                else:                    
                    genename = linedata[0].upper()
                    if genename.split("_")[0] == "CRID":
                        genename = "COMMUNICATIONROUEID_"+ genename.split("_")[1]
                    # print(genename)
                    if genename not in genelist:
                        genelist.append(genename)
                    del(linedata[0])
                    for i in range(len(samplelist)):
                        allnamesamples[samplelist[i]][genename] = float(linedata[i])
    return allnamesamples,genelist

def loadfilewithoutname(inputfilepath, transpose=False):
    
    '''
    Load file from txt file

    First row should be sample name
    First column should be gene name
    other wise transpose should be True
    input 
        inputfilepath: Path to file
        transpose: Does this file need transpose. Default is transpose
    return a 2-D list. [number][genename] = float value 
    '''

    allsamples=[]
    with open(inputfilepath,'r') as allsamplefile:
        istitle=True
        if transpose:
            genelist = [] 
            for line in allsamplefile.readlines():
                linedata = line.strip().split("\t")
                if istitle:
                    del(linedata[0])
                    for genename in linedata:
                        genelist.append(genename.upper())                        
                    istitle=False
                else:
                    sampledetails ={} 
                    del(linedata[0])
                    for i in range(len(genelist)):                        
                        sampledetails[genelist[i]] = (float(linedata[i]))
                    allsamples.append(sampledetails)                    
        else:
            for line in allsamplefile.readlines():
                linedata = line.strip().split("\t")
                if istitle:
                    del(linedata[0])
                    for sample in linedata:
                        allsamples.append({})
                    istitle=False
                else:                    
                    genename = linedata[0].upper()
                    del(linedata[0])
                    for i in range(len(linedata)):
                        allsamples[i][genename] = (float(linedata[i]))
    return allsamples

def  seperatesampleslist(numberofgroup,allsamples,israndom=True):
    '''
    This function is to split data list randomly into subgroup

    input 
        numberofgroup: Path to file
        allsamples: Alist of sample  [number][genename] = float value 
        israndom: randomly split or not. Default is True
    
    return  a list [numberofgroup][number][genename]

    '''

    shuffledname=list(allsamples.keys())

    if israndom:
        random.shuffle(shuffledname)

    seperatedsamples=[]

    for i in range(numberofgroup):
        seperatedsamples.append({})

    index =-1
    for samplename in shuffledname:
        if index==numberofgroup-1:
            index =0
        else:
            index+=1
            
        seperatedsamples[index][samplename]=allsamples[samplename]

    return seperatedsamples

def getGeneVector(allsamples, genename):

    geneVector = []
    
    for samplename, sample in allsamples.items():
        if genename in sample.keys():
            geneVector.append(float(sample[genename]))
        else:
            geneVector.append('Nan')

    return geneVector

def getsampleVectors(sample, genenamelist):

    sampleVector = []
    
    for gene in genenamelist:
        if gene in sample.keys():
            sampleVector.append(float(sample[gene]))
        else:
            sampleVector.append(float('Nan'))

    
    return sampleVector

def writesampletofile(allsamples,outputfile):

    '''
    Write samples list into file

    input 
        allsamples: Alist of sample  [number][genename] = float value 
        outputfile: file path
    
    return  NA

    '''

    allgene = allsamples[0].keys()

    with open(outputfile,'w') as output:
        title = ""
        for genename in allgene:
            title+=genename.upper()+"\t"
        title = title[:-1]+"\n"
        output.write(title)

        for sample in allsamples:
            pline = ""
            for genename in allgene:
                pline+=sample[genename.upper()]+"\t"
            pline = pline[:-1]+"\n"
            output.write(pline)

def writenamedsampletofile(allsamples,outputfile):

    '''
    Write samples list into file

    input 
        allsamples: Alist of sample  [name][genename] = float value 
        outputfile: file path
    
    return  NA

    '''
    if len(allsamples)>0:
        samplelist = list(allsamples.keys())
        genelist = list(allsamples[samplelist[0]].keys())

        with open(outputfile,'w') as output:
            title = "GENE\t"
            for samplename in samplelist:
                title+=samplename.upper()+"\t"
            title = title[:-1]+"\n"
            output.write(title)

            for genename in genelist:
                pline = genename.upper()+"\t"
                for samplename in samplelist:
                    pline+=str(allsamples[samplename][genename])+"\t"
                pline = pline[:-1]+"\n"
                output.write(pline)

def computesumofsample(sample):
    return sum(sample.values())

def changepotionofsample(sample,weightfactor):
    ratiosample = {}
    sumofsample = computesumofsample(sample)
    for genename,value in sample.items():
        ratiosample[genename] = value/sumofsample*weightfactor
    return ratiosample

def changepotionofsamples(samples,weightfactor=10000):

    newsamples={}
    for samplename,sample in samples.items():
        newsamples[samplename] = changepotionofsample(sample,weightfactor)
    return newsamples

def getthrankedsamplelistofsum(samples):
    sumdict = {}
    for samplename,sample in samples.items():
        sumdict[samplename]= computesumofsample(sample)
    return  sumdict

def pickthethresholdofnamedsamples(sumdict,bottompercent,toppercent):

    rankedsumdict = {k: v for k, v in sorted(sumdict.items(), key=lambda item: item[1],reverse=True)}

    numberofsample = len(sumdict)
    selectedsamplenames = []
    numberofbottom = numberofsample*bottompercent
    numberoftop = numberofsample*toppercent

    index = 1
    for samplename in list(rankedsumdict.keys()):
        
        if index>numberofbottom and index < numberoftop:
            selectedsamplenames.append(samplename)
        index +=1
    return selectedsamplenames

def getsubgroupfromsamples(namedsamples,selectednamelist):
    selectedgroups= {}
    
    for samplename in selectednamelist:
        if samplename in namedsamples.keys():
            selectedgroups[samplename] = namedsamples[samplename]
    
    return selectedgroups

def computesubgroupsratio(samples):
    
    ratiodict = {}
    sumdict = {}
    numberofsamples = len(samples)
    totalsum = 0
    for samplename,sample in samples.items():
        for genename, value in sample.items():
            if genename not in sumdict.keys():
                sumdict[genename]=0
            sumdict[genename]+= value
        totalsum+=computesumofsample(sample)
    
    for genename,sumvalue in sumdict.items():
        ratiodict[genename] = sumvalue/totalsum
    
    return ratiodict

def changeferencesample(sample, referenceratio):
    ratiosample = changepotionofsample(sample,1)
    newsample = {}
    '''
    if test  == 0 
       give 0 
    if reference == 0
       give 5
    '''
    for genename,value  in ratiosample.items():
        if ratiosample[genename]==0:
            if referenceratio[genename]==0:
                newsample[genename] = 0
            else:
                newsample[genename] = -15
        else:
            if referenceratio[genename]==0:
                newsample[genename] = 15
            else:
                newsample[genename]= math.log2(ratiosample[genename]/referenceratio[genename])
    
    return newsample

def referencechangesamples(originalnamedsamples,bottompercent,toppercent):

    sumdict = getthrankedsamplelistofsum(originalnamedsamples)
    selectedsamplenames = pickthethresholdofnamedsamples(sumdict,bottompercent,toppercent)
    selectedsamples = getsubgroupfromsamples(originalnamedsamples,selectedsamplenames)

    referencedict = computesubgroupsratio(selectedsamples)

    # here first we try to take log 2 of ratio between test and reference
    newsamples = {}
    for samplename,sample in originalnamedsamples.items():
        newsamples[samplename] = changeferencesample(sample, referencedict)
    
    return newsamples

def samplestomatrix(samples):

    matrix = []
    for samplename,sample in samples.items():
        genevalues = list(sample.values())
        matrix.append(genevalues)
    return matrix

def computethecorrelationofselectedgenes(samples,selectedgenes):
    selecteddict = {}
    correlationdict={}


    for genename in selectedgenes:
        selecteddict[genename] =  getGeneVector(samples, genename)
    
    
    corcombinations = list(combinations(selectedgenes, 2))

    for corcombination in corcombinations:
        genelist1=selecteddict[corcombination[0]]
        genelist2=selecteddict[corcombination[1]]

        corr, _ = pearsonr(genelist1 ,genelist2 )

        labelname = corcombination[0]+','+corcombination[1]
        # print(labelname)
        correlationdict[labelname]=corr
    
    return correlationdict

def writeCorrelationtoFile(correlationdict,outputfile):

    correlationdict={k: v for k, v in sorted(correlationdict.items(), key=lambda item: item[1],reverse=True)}
    with open(outputfile,'w') as output:
        for label,correlstionvalue in correlationdict.items():
            
            output.write(label.replace(',','\t')+'\t'+str(correlstionvalue)+'\n')
    
def gettherankofgenesum(samples,genelist):
    rankedgene = {}
    for genename in genelist:
        rankedgene[genename]=sum(getGeneVector(samples,genename))
    
    rankedgene={k: v for k, v in sorted(rankedgene.items(), key=lambda item: item[1],reverse=True)}
    
    
    return rankedgene 

def removeselectedgenefromsample(sample,selectedgenelist):
    newsample = {}
    for genename, genevalue in sample.items():
        if genename not in selectedgenelist:
            newsample[genename] = genevalue
    
    return newsample

def removeselectedgenefromsamples(samples,selectedgenelist):
    newsamples = {}
    for samplename,sample in samples.items():
        newsamples[samplename] = removeselectedgenefromsample(sample,selectedgenelist)
      
    return newsamples 

def selectedgenefile(selectedfile):
    selectedgenelist = []
    with open(selectedfile,'r') as selectedlist:
        for line in selectedlist.readlines():
            selectedgenelist.append(line.strip())
    
    return selectedgenelist

def selectedgenesbycorrelation(correlationdict,numberofselected,prefix = "RP"):

    correlationdict={k: v for k, v in sorted(correlationdict.items(), key=lambda item: item[1],reverse=True)}

    selectedgenenlist = []

    for keypairs in correlationdict.keys():
        [gene1,gene2]=keypairs.split(',')
        if gene1.startswith(prefix) and gene1 not in selectedgenenlist and len(selectedgenenlist)<=numberofselected:
            selectedgenenlist.append(gene1)
        if gene2.startswith(prefix) and gene2 not in selectedgenenlist and len(selectedgenenlist)<=numberofselected:
            selectedgenenlist.append(gene2)
    return selectedgenenlist

def loadcorrelationfile(inputfile):
    correlationdict  = {}

    with open(inputfile,'r') as input:
        for line in input.readlines():
            linedata = line.strip().split('\t')
            label = linedata[0]+','+linedata[1]
            value = float(linedata[2])
            correlationdict[label]=value
    
    return correlationdict

def writeselectedgeneintofile(filepath,genelist):
    with open(filepath,'w') as outfile:
        for gene in genelist:
            outfile.write(gene+'\n')

def findCommonlistinlist(genelistlist):

    result = set(genelistlist[0])
    for s in genelistlist[1:]:
        result.intersection_update(s)
    return result

def findCommongeneintwolist(genelistlist1,genelist2):

    return  set(genelistlist1). intersection(set(genelist2))

def getelectedgenefromsample(sample,selectedgenelist):
    newsample = {}
    for genename, genevalue in sample.items():
        if genename in selectedgenelist:
            newsample[genename] = genevalue
    
    return newsample

def getselectedgenefromsamples(samples,selectedgenelist):
    newsamples = {}
    for samplename,sample in samples.items():
        newsamples[samplename] = getelectedgenefromsample(sample,selectedgenelist)
      
    return newsamples 

def loadlabelfile(filepath,sampenamecol,labelnamecol,sep='\t'):

    labeldict={}
    istitle=True
    with open(filepath,'r') as inputfile:
        for line in inputfile.readlines():
            if istitle:
                istitle=False
            else:
                linedata = line.strip().split(sep)
                samplename = linedata[sampenamecol].upper()
                labelname = linedata[labelnamecol].upper().replace("_","-")
                
                if labelname not in labeldict.keys():
                    labeldict[labelname] = []
                labeldict[labelname].append(samplename)
    
    return  labeldict

def writePostanalysisVariant(Variants,outputpath):
    with open(outputpath,'w') as outputfile:
        for gene1,gene2list in Variants.items():
            for gene2 in gene2list:
                outputfile.write(gene1+'\t'+gene2+'\n')

def writefilteredweighttofile(Weights,outputfile):
    '''
    write subggroup ranked weight to output file
    inputs 
        totalweight[targetname][casualgene]=[rank1, rank2, rank3, rank4...]
    outputs
        a file that record the weight
    ''' 
    with open(outputfile, 'w') as output:
        for target, causaltagetdetails in Weights.items():
            for casualgene,corvalue in causaltagetdetails.items():
                printline = target+'\t'+casualgene +'\t'+str(corvalue)                
                output.write(printline+'\n')

def writeCorrelationdicttofile(outputfilepath, correlationdict):

    with open(outputfilepath, 'w') as output:
        for label,corr in correlationdict.items():
            printline = label+'\t'+str(corr)                
            output.write(printline+'\n')

def computeavgerofsamples(samples,genelist):
    averagedict={}
    for genename in genelist:
        genevector = getGeneVector(samples,genename)
        averagevalue = statistics.mean(genevector)
        averagedict[genename] = averagevalue

    return averagedict

def computelogofsample(sample,controlsample):
    log2dict = {}
    for genename,value in sample.items():
        controlvalue = controlsample[genename]
        if value >0 :
            if controlvalue>0:
                log2dict[genename] = math.log2(value/controlvalue)
            else:
                log2dict[genename] = 0

        else:
            if controlvalue>0:
                log2dict[genename] = 0
            else:
                log2dict[genename] = 0

    return log2dict

def computelogofsamples(samples,controlsample):
    log2samples = {}
    for samplename,sample in samples.items():
        log2samples[samplename]=computelogofsample(sample,controlsample)
    
    return log2samples



def calculategenewithsumcorrelation(samples,genelist):



    sumcorrdict= {}

    sumlist = []
    for samplename, sample in samples.items():
        sumlist.append(computesumofsample(sample))
    
    for gene in genelist:
        genevector = getGeneVector(samples,gene)
        corr, _ = pearsonr(genevector ,sumlist )
        sumcorrdict[gene] = corr
    
    return sumcorrdict


def selectedgenesbycorrelationsum(correlationdict,numberofselected,prefix = "RP"):

    correlationdict={k: v for k, v in sorted(correlationdict.items(), key=lambda item: item[1],reverse=True)}

    selectedgenenlist = []

    for gene in correlationdict.keys():
        if gene.startswith(prefix) and gene not in selectedgenenlist and len(selectedgenenlist)<=numberofselected:
            selectedgenenlist.append(gene)
        
    return selectedgenenlist


def nthrootofm(a,n):
    return pow(a,(1/n))


def writingwithrow(samples,outputfile):
    genelist = samples[list(samples.keys())[0]].keys()
    with open(outputfile,'w') as output:
        headerline = "Cellname\t"+'\t'.join(genelist)+'\n'
        output.write(headerline)
        for cellname, cellvalue in samples.items():
            line = '\t'+'\t'.join(list(cellvalue.values()))+'\n'
            output.write(line)




def generaterandomdataset(rawfilepath,outputfile):
    randomizdataset = {}
    celldataset,genelist = loadfilewithname(rawfilepath,transpose=True)
    
    for cellname in list(celldataset.keys()):
        print("random",cellname)
        newsample = {}
        selectedvaluelist = list(random.choice(list(celldataset.values())).values())
        for genename in genelist:
            newsample[genename] = str(random.choice(selectedvaluelist))
        randomizdataset[cellname]= newsample
    
        randomizdataset[cellname]= newsample
    writingwithrow(randomizdataset,outputfile)



def getrootvalue(samples,genenameindex,lappa=0.0001):
    routesubunit = 0
    totalvalue  = 1

    
    for samplename,sample in samples.items():
        genevalue =float(sample[genenameindex])+lappa
        if genevalue>0:
            totalvalue*=genevalue
            routesubunit+=1

    if routesubunit>0:
        return nthrootofm(totalvalue,routesubunit)
    else:
        return 0




    




















