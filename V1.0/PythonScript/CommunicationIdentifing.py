
import GeneExpressionOperate

from itertools import combinations



def CalculateGroupRootValue(samples, routelist):
    groupsamples={}
    
    for i  in range(len(routelist)) :
        CommunicationRouteRootValue = GeneExpressionOperate.getrootvalue(samples,i,lappa=1)
        groupsamples[routelist[i]] = CommunicationRouteRootValue
    
    return groupsamples
        
        

def CalculatePobablility(CellgroupSamples,routelist):

    Cellgroupsnamelist = list(CellgroupSamples.keys())
    totallist = []
    # combinations =  combinations(Cellgroupsnamelist, 2)
    Pdict = {}
    for i in range(len(routelist)):
        
        CLRroutename = routelist[i]
        if "CLR" in CLRroutename:
            j = i +1
            CRRroutename = routelist[j]
            routename = CLRroutename.split("-")[0]
            print(routename)
            for Lgroupname in Cellgroupsnamelist:
                # LARA = CellgroupSamples[Lgroupname][CLRroutename] * CellgroupSamples[Lgroupname][CRRroutename]
                for Rgroupname in Cellgroupsnamelist:
                    # LARBlist=[]
                    if Lgroupname+'_'+Rgroupname not in Pdict.keys():
                        Pdict[Lgroupname+'_'+Rgroupname]={}
                    
                    LARB = CellgroupSamples[Lgroupname][CLRroutename] * CellgroupSamples[Rgroupname][CRRroutename]
                    
                    Pdict[Lgroupname+'_'+Rgroupname][routename] = LARB
                    totallist.append(LARB)
                    # if not Rgroupname == Lgroupname:

                    #     LARB = CellgroupSamples[Lgroupname][CLRroutename] * CellgroupSamples[Rgroupname][CRRroutename]
                    #     P = (LARB/100)*(1-LARA/100)
                    #     Pdict[Lgroupname+'_'+Rgroupname][routename] = P
                        
                    # # else:
                    #     P = (LARA/100)*(LARA/100)
                    #     Pdict[Lgroupname+'_'+Rgroupname][routename] = P
    maxvalue = max(totallist)
    newPdict={}
    for samplename,routenamevalues in Pdict.items():       
        newPdict[samplename]={}
        for routename,routevalue in routenamevalues.items():
            print(samplename,routename)            
            newPdict[samplename][routename]= str(round(routevalue/maxvalue,2))
       
    return newPdict


def CalculatePobablility_v2(CellgroupSamples,routelist):

    Cellgroupsnamelist = list(CellgroupSamples.keys())
    totallist = []
    # combinations =  combinations(Cellgroupsnamelist, 2)
    Pdict = {}
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
                    if Lgroupname+'_'+Rgroupname not in Pdict.keys():
                        Pdict[Lgroupnamenew+'_'+Rgroupnamenew]={}
                    
                    LARB = CellgroupSamples[Lgroupname][CLRroutename] * CellgroupSamples[Rgroupname][CRRroutename]
                    LBRB = CellgroupSamples[Rgroupname][CLRroutename] * CellgroupSamples[Rgroupname][CRRroutename]

                    Pdict[Lgroupnamenew+'_'+Rgroupnamenew][routename] = (LARB+1)/(LARA+LBRB+1) *10
                    # totallist.append(LARB)
                    # if not Rgroupname == Lgroupname:

                    #     LARB = CellgroupSamples[Lgroupname][CLRroutename] * CellgroupSamples[Rgroupname][CRRroutename]
                    #     P = (LARB/100)*(1-LARA/100)
                    #     Pdict[Lgroupname+'_'+Rgroupname][routename] = P
                        
                    # # else:
                    #     P = (LARA/100)*(LARA/100)
                    #     Pdict[Lgroupname+'_'+Rgroupname][routename] = P
    # maxvalue = max(totallist)
    # newPdict={}
    # for samplename,routenamevalues in Pdict.items():       
    #     newPdict[samplename]={}
    #     for routename,routevalue in routenamevalues.items():
    #         print(samplename,routename)            
    #         newPdict[samplename][routename]= str(round(routevalue/maxvalue,2))
       
    return Pdict               
    

def calculateCommunicationPability(CommunicationRouteScoreSamples, sampleslabel,routelist,subtmpfile):
    CellgroupSamples = {}

    for celltypename, cellsamplenamelist in sampleslabel.items():
        subsamples = GeneExpressionOperate.getsubgroupfromsamples(CommunicationRouteScoreSamples,cellsamplenamelist)
        CommunicationsubsamplesRootVectors= CalculateGroupRootValue(subsamples, routelist)
        CellgroupSamples[celltypename] = CommunicationsubsamplesRootVectors
        print(celltypename)

    GeneExpressionOperate.writenamedsampletofile(CellgroupSamples,subtmpfile)
    
    return CalculatePobablility(CellgroupSamples,routelist)



def calculateCommunicationPability_v2(CommunicationRouteScoreSamples, sampleslabel,routelist,subtmpfile,celltypefolder):
    CellgroupSamples = {}

    for celltypename, cellsamplenamelist in sampleslabel.items():
        subsamples = GeneExpressionOperate.getsubgroupfromsamples(CommunicationRouteScoreSamples,cellsamplenamelist)

        with open(celltypefolder+celltypename+'_communication.txt','w') as outputfile:
            outputfile.write("Cellname\t"+"\t".join(routelist)+'\n')
            for samplename,samplevalue in subsamples.items():
                strsamplevalue =[]
                for a in samplevalue:
                    strsamplevalue.append(str(a))
                outputfile.write(samplename+"\t"+"\t".join(strsamplevalue)+'\n')

        # GeneExpressionOperate.writenamedsampletofile(subsamples,celltypefolder+celltypename+'_communication.txt')

        CommunicationsubsamplesRootVectors= CalculateGroupRootValue(subsamples, routelist)
        CellgroupSamples[celltypename] = CommunicationsubsamplesRootVectors
        print(celltypename)
    GeneExpressionOperate.writenamedsampletofile(CellgroupSamples,subtmpfile)
    return CalculatePobablility_v2(CellgroupSamples,routelist)





if __name__=="__main__":




    # outputfolder= './Dataset/Bone/Mouse bone study1/Raw_processed/'
    
     
    # # celllogfilepath = outputfolder+"/preprocessed_log.txt"

    # Communication_bonesamplefile = outputfolder+"/Communication_preprocessed_log_v2.txt"
    # outputfolder= './Dataset/Bone/Mouse bone study1/Raw_processed/'

    # # outputfolder= './Dataset/Bone/GSE152285/RAW/'
    # # # celllogfilepath = outputfolder+"/GSE152285_preprocessed_log.txt"
     
    # # Communication_bonesamplefile = './Dataset/Bone/GSE152285/RAW/Communication_preprocessed_log_testcorrect.txt'

    # tempfile = outputfolder+"Communication_preprocessed_log_testcorrect_tmp.txt"

    # # outputfolder= './Dataset/Bone/GSE152285/RAW/'

    

    
    # # celllogfilepath =outputfolder+"/GSE152285_preprocessed_log.txt"
     
    # Communication_bonesamplefile = outputfolder+"/Communication_preprocessed_log_v2.txt"


    
    outputfolder= './Dataset/Bone/GSE152285/RAW_test6/'
    celllogfilepath = outputfolder+"/preprocessed.txt"
    tempfile = outputfolder+"Communication_preprocessed_log_testcorrect_tmp.txt" 
    Communication_bonesamplefile = outputfolder+'/Communication_preprocessed.txt'

    # Communication_bonesamplefile = outputfolder+"/Communication_probability_preprocessed_log_v3.txt"

    # load cell label
    # labelfilepath = "./Dataset/Bone/Mouse bone study1/mc_compliant_metadata_example.tsv"
    # sampenamecol = 0
    # labelnamecol = -2


    labelfilepath = "./Dataset/Bone/GSE152285/meta_annotation_revised.txt"
    sampenamecol = 0
    labelnamecol = -2



    cellleidenlabel = GeneExpressionOperate.loadlabelfile(labelfilepath,sampenamecol,labelnamecol,"\t")

    print("Done load label")
    # load communication scoring 
    Communicationfilepath =  Communication_bonesamplefile

    communication_routescoredict = {}

    with open(Communicationfilepath,'r') as Celldataset:
        # load first row as genename
        
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
                
    # CommuncationCells,routelist = GeneExpressionOperate.loadfilewithname(Communicationfilepath,transpose=True)
    print("Done load communication")

    # pobabilityoutputfile = outputfolder+"/Communication_probability_GSE152285_preprocessed_log_v3.txt"

    pobabilityoutputfile = outputfolder+"/Communication_probability_preprocessed_log_v5.txt"



    # Pobalityresult = calculateCommunicationPability(communication_routescoredict, cellleidenlabel,routelist)
    Pobalityresult = calculateCommunicationPability_v2(communication_routescoredict, cellleidenlabel,routelist,tempfile,outputfolder)

    print("Done load Calculating")

    GeneExpressionOperate.writenamedsampletofile(Pobalityresult,pobabilityoutputfile)