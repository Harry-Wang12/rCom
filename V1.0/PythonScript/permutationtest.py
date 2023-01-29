
import random
import CommunicationIdentifing
import GeneExpressionOperate


def CalculateGroupRootValue(samples, routelist):
    groupsamples={}
    
    for i  in range(len(routelist)) :
        CommunicationRouteRootValue = GeneExpressionOperate.getrootvalue(samples,i,lappa=1)
        groupsamples[routelist[i]] = CommunicationRouteRootValue
    
    return groupsamples
        

def CalculatePobablility_per(CellgroupSamples,routelist):

    Cellgroupsnamelist = list(CellgroupSamples.keys())
    totallist = []
    # combinations =  combinations(Cellgroupsnamelist, 2)
    Pdict = []
    for i in range(len(routelist)):
        
        CLRroutename = routelist[i]
        if "CLR" in CLRroutename:
            j = i +1
            CRRroutename = routelist[j]
            routename = CLRroutename.split("-")[0]
            print(routename)
            for Lgroupname in Cellgroupsnamelist:
                Lgroupnamenew = Lgroupname.replace("_","-")
                LARA = CellgroupSamples[Lgroupname][CLRroutename] * CellgroupSamples[Lgroupname][CRRroutename]
                for Rgroupname in Cellgroupsnamelist:
                    Rgroupnamenew = Rgroupname.replace("_","-")
                    # LARBlist=[]
                   
                    
                    LARB = CellgroupSamples[Lgroupname][CLRroutename] * CellgroupSamples[Rgroupname][CRRroutename]
                    LBRB = CellgroupSamples[Rgroupname][CLRroutename] * CellgroupSamples[Rgroupname][CRRroutename]

                    Pdict.append((LARB+1)/(LARA+LBRB+1) *10)
      
    return Pdict

def calculateCommunicationPability_per(CommunicationRouteScoreSamples, sampleslabel,routelist,subtmpfile,celltypefolder):
    CellgroupSamples = {}

    for celltypename, cellsamplenamelist in sampleslabel.items():
        subsamples = GeneExpressionOperate.getsubgroupfromsamples(CommunicationRouteScoreSamples,cellsamplenamelist)
        CommunicationsubsamplesRootVectors= CalculateGroupRootValue(subsamples, routelist)
        CellgroupSamples[celltypename] = CommunicationsubsamplesRootVectors
        print(celltypename)
    GeneExpressionOperate.writenamedsampletofile(CellgroupSamples,subtmpfile)
    return CalculatePobablility_per(CellgroupSamples,routelist)




if __name__=="__main__":



    outputfolder= './Dataset/Bone/GSE152285/RAW_test6/'
    Communication_bonesamplefile = outputfolder+"/Communication_preprocessed.txt"

    permutationtest = outputfolder+"/Communication_probability_preprocessed_log_randompermutationtest.txt"

    tempfile = outputfolder+"/Communication_probability_preprocessed_log_randompermutationtest_tmp.txt"

    # cellsupgroup1 = {}
    # cellsupgroup2 = {}

    testnumber = 10
    k = 100

    Communicationfilepath =  Communication_bonesamplefile

    communication_routescoredict = {}

    with open(Communicationfilepath,'r') as Celldataset:
        istitle = True
        
        needgeneindex = []

        isstop = False

        # genedict = manager.dict()
        while True:
            line = Celldataset.readline()
            if line:
                if istitle:
                    
                    genelist = line.strip().split('\t')[1:]
                    routelist = []
                                    
                    for i in range(len(genelist)):
                        routename  = genelist[i]
                        routelist.append(routename)

                    istitle = False
                    i=0
                else:
                    i+=1
                    print(i)
                    samplevaluelist = []
                    linedata = line.strip().split('\t')
                    samplename = linedata[0]
                    
                    del(linedata[0])
                    for value in linedata:
                        samplevaluelist.append(value)
                    communication_routescoredict[samplename] = samplevaluelist
            else:
                break
    
    samplenamelist = list(communication_routescoredict.keys())
    labellist = {}
    for i in range(testnumber):
        labelname = str(i)
        labellist[labelname]= random.choices(samplenamelist, k = 100)

    Pobalityresultlist = calculateCommunicationPability_per(communication_routescoredict, labellist,routelist,tempfile,outputfolder)


    print("Done load Calculating")

    # GeneExpressionOperate.writenamedsampletofile(Pobalityresult,permutationtest)

   



    summaryfile = outputfolder+"summary.txt"
    newsummaryfile = outputfolder+"summary_pvalue.txt"

    with open(summaryfile,'r') as summary:
        with open(newsummaryfile,'w') as summaryout:
            istitle = True
            for line in summary.readlines():
                if istitle:
                    summaryout.write(line.strip()+'\tpval\n')
                    istitle = False
                else:
                    linedata = line.strip().split('\t')
                    communval = float(linedata[3])
                    pvalue = 0
                    for value in Pobalityresultlist:
                        if value <= communval:
                            pvalue+=1
                    pvalue =1- pvalue/len(Pobalityresultlist)

                    summaryout.write(line.strip()+'\t'+str(pvalue)+'\n')