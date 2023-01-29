





if __name__ == "__main__":


    getsubgroup = False

    outputfolder= './Dataset/Bone/GSE152285/RAW_test6/'
    samplelist=[]
    with open(outputfolder+"/obs.csv",'r') as samplenamefile:
        title = True
        for line in samplenamefile.readlines():
            if title:
                title=False
            else:
                samplename = line.strip().split(',')[0]

                samplelist.append(samplename)

    
    celltypedsamples = []

    labelfilepath = "./Dataset/Bone/GSE152285/meta_annotation_revised.txt"
    sampenamecol = 0
    labelnamecol = 1

    with open(labelfilepath,'r') as labelfile:
        for line in labelfile.readlines():
            samplename = line.strip().split('\t')[0]
            celltypedsamples.append(samplename)

    lst3 = [value for value in samplelist if value in celltypedsamples]

    print(len(lst3))

    if getsubgroup:
        intputfile = outputfolder+'preprocessed.txt'
        outputfile = outputfolder+'sub_preprocessed.txt'
        title = True
        with open(intputfile,"r") as input:
            with open(outputfile,'w') as output:
                for line in  input.readlines():
                    if title:
                        output.write(line)
                        title = False
                    else:
                        linedata = line.strip().split('\t')
                        if linedata in celltypedsamples:
                            output.write(line)

