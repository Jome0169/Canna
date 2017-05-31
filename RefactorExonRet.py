import os
import sys
import getopt
import itertools
import math
from GeneClass import GeneObj
from datetime import datetime
from itertools import combinations


def GenomeReader(GenomeFile):
    """
    Arg: Takes in Genome File
    Rtrns: Returns a dictionary, Genome Scaffolds. 
    
    Keys genomic scaffold names being the
    keys - and the actual sequence being the value. 
    """
    GenomeScaffolds = {}
    with open(GenomeFile, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                NamedSeq = line.replace('>', '')
                GenomeScaffolds[NamedSeq] = ""
            else:
                GenomeScaffolds[NamedSeq] += line
        return GenomeScaffolds


def BlastFilereader(BLastfile):
    # BlastFilereader will return a list of list from a blast file with the
    #above commented 7 header feature from a BLAST table."""
    HitsNStuff = []
    with open(BLastfile, "r") as f:
        for line in f:
           HitsNStuff.append(line.strip("\n").split('\t'))
    return HitsNStuff

def Scafs2Dict(ListofBlast):
    #Creates a dictionary with all uniq scaffolds name from the assembly being
    #Keys in the dictionary
    
    #Below is an example of stuff
    #002261F|quiver
    #{}

    UniqScaffolds = {}
    for item in ListofBlast:
        if item[1] not in UniqScaffolds:
            UniqScaffolds[item[1]] = {}
    return UniqScaffolds



def BlastParser(ListofBlast, DictofUniqScaffolds):
    #Takes in the BLAST file again, and then appends these items from the blast
    #table to the dictionary based off of the Scaffolds that these inidivuals
    #hit to. 
    for item in ListofBlast: #Creates nested dictiomary item Fasta header Seq K
        DictofUniqScaffolds[item[1]][item[0]] = []
    for item in ListofBlast: # Adds Blast hits to 
        DictofUniqScaffolds[item[1]][item[0]].append(item)
    return DictofUniqScaffolds


def DictEditor(DictGeneObj):
    """TODO: Docstring for DictEditor.

    :DictGeneObj: TODO
    :returns: TODO

    """

    def ScaffoldHitRange(ListofList):
        """TODO: Docstring for ScaffoldHitRange.
        
        :ListofList: TODO
        :returns: TODO
        
        """

        DicrectionInfo = CheckDirection(ListofList)

        if DicrectionInfo == 0:
            LargestScafHitRegion = max(ListofList, key=lambda x: int(x[5]))
            SmallestScafHitRegion = min(ListofList, key=lambda x: int(x[4]))
            SmallestInt = int(SmallestScafHitRegion[4])
            LargestInt = int(LargestScafHitRegion[5])
            
            if SmallestInt < LargestInt:
                ScaffoldRange = (SmallestInt,
                        LargestInt)
                return ScaffoldRange
            
        elif DicrectionInfo == 1:
            LargestScafHitRegion = max(ListofList, key=lambda x: int(x[4]))
            SmallestScafHitRegion = min(ListofList, key=lambda x: int(x[5]))
            SmallestInt = int(SmallestScafHitRegion[5])
            LargestInt = int(LargestScafHitRegion[4])

            if SmallestInt < LargestInt:
                ScaffoldRange = (SmallestInt,
                        LargestInt)
                return ScaffoldRange

        elif DicrectionInfo == 2:
            print "ERROROOROROR"
            for item in ListofList:
                print item

    def CheckDirection(ListofList):
        """TODO: Docstring for .

        :ListofList: TODO
        :returns: TODO

        """
        Posotive = 0
        Negative = 0
        
        for blastresult in ListofList:
            if int(blastresult[4]) < int(blastresult[5]):
                Posotive += 1
            elif int(blastresult[4]) > int(blastresult[5]):
                Negative += 1

        if Posotive > 0:
            return 0
        elif Negative > 0:
            return 1
        elif Posotive > 0 and Negative > 0:
            return 2



    
    def OverlapCheckandElimination(ListofList1, TupleLength1, ListofList2, \
            TupleLength2):
        """TODO: Docstring for OverlapCheckandElimination.

        :tuplepair: TODO
        :returns: TODO

        """

        count1 = len(ListofList1)
        count2 = len(ListofList2)
        

        #Maybe come back and calculate AVG exon distance instead of "length"
        Length1 = int(TupleLength1[1]) - int(TupleLength1[0])
        Length2 = int(TupleLength2[1]) - int(TupleLength2[0])

        EScoreCalc1 = 0
        EScoreCalc2 = 0

        for item in ListofList1:
            if float(item[7]) == 0:
                EScoreCalc1 += 0
            else:
                EScoreCalc1 += math.log(float(item[7]))
        
        EScoreCalc1 = EScoreCalc1 / count1
        
        for item in ListofList2:
            if float(item[7]) == 0:
                EScoreCalc2 += 0
            else:
                EScoreCalc2 += math.log(float(item[7]))
        
        EScoreCalc2 = EScoreCalc2 / count2


        #LOGIC TO DELET SHIT
        #Length = Lenght of sequence
        #count =  Number of Exons
        #Escore  = sum of Escores

        #The below logical returns the items that are of low quality and will
        #be removed in later processing steps
        if count1 > count2:
            return ListofList2[0][0]

        elif count2 > count1:
            return ListofList1[0][0]

        elif count1 == count2 and float(EScoreCalc1) < float(EScoreCalc2):
            return ListofList2[0][0]

        elif count1 == count2 and float(EScoreCalc2) < float(EScoreCalc1):
            return ListofList1[0][0]

        elif count1 == count2 and EScoreCalc1 == EScoreCalc1 and Length1 == Length2:
            return ListofList1[0][0]
        
        elif count1 == count2 and EScoreCalc1 == EScoreCalc2 and Length1 > Length2:
            return ListofList2[0][0]

        elif count1 == count2 and EScoreCalc1 == EScoreCalc2 and Length1 < Length2:
            return ListofList1[0][0]

        else: 
            print "ERROR"
           



    for Scaffold, Gene in DictGeneObj.iteritems():
        lengths = [len(v) for v in Gene.values()] #will be useful
    
        if len(lengths) == 1:
            pass
        
        else:
            GeneToBeRemoved = []
            KeyList = Gene.keys()
            PossiblePairs = list(itertools.combinations(KeyList, 2))
            
            for pairing in PossiblePairs:
                if pairing[0] in Gene.keys() and pairing[1] in Gene.keys():
                    CheckRange1 = ScaffoldHitRange(Gene[pairing[0]])
                    CheckRange2 = ScaffoldHitRange(Gene[pairing[1]])
                    #Below checks Overlap

                    if CheckRange1[0] <= CheckRange2[1] and CheckRange2[0] <= \
                    CheckRange1[1]:
                        RemoveThis = OverlapCheckandElimination(Gene[pairing[0]], CheckRange1, \
                                Gene[pairing[1]], CheckRange2)
                        GeneToBeRemoved.append(RemoveThis)
                    else:
                        pass
                        
    
            Unique = set(GeneToBeRemoved)
            for LowQualGene in Unique:
                if LowQualGene != None:
                    del Gene[LowQualGene]
                else:
                    pass


    return DictGeneObj
    


def PullingOutLowExonGenes(arg1):
    """TODO: Docstring for PullingOutLowExonGenes.

    :arg1: Takes in list of blast hits that have less than 5 exons. Individuals
    with seq less than 400 Nts in length are thrown out. None type is returned
    in that case.  The idea here is that they might be psuedo genes or just have a few segments.
    However there are quite a few indivuduals with large hit boxes. So these
    will be pulled

    """


    if len(arg1) == 1 and (int(arg1[0][3]) - int(arg1[0][2])) <= 100:
        #Length Too Small
        return 1
    elif len(arg1) == 1 and (int(arg1[0][3]) - int(arg1[0][2])) >= 400:
        #Probably Good Gene. return arg
        return 2
    elif len(arg1) > 1: # If greater than length 1, begin more filtering 
        y = Filter2(arg1)
        return y

def Filter2(arg2):
    """ If there is more than 1 exon we have to spend some time
    sorting the exons as well as looking for duplication. Takes in list,
    ensure total sum is greater than 800. Checks for overlap between the
    genes using protein query hit lovation (if both start with one)

    :arg2: List of list. Genes with less than 5 exons but more than 1. If
    total len is less than 400 Seq is deleted and Nonetype retruned
    :returns: None type or a list of list

    """
    NucLenSum = 0
    for item in arg2: # Calculate how much Nucleotide seq there is 
        SeqLen = (int(item[3]) - int(item[2]))
        NucLenSum += SeqLen 
    if NucLenSum >= 400: # if Greater than 800 nuceltodies (260 AA),proceed
        return 2
    else:
        #Gene too short. Delete Argument. Will return None
        return 1


def FindDuplications(arg1):
    """FindDuplications takes in genomic hits and searches for a
    gene/transcript hitting to a genomic scaffold multiple times. If this is
    the case we then analyze weather the Start of the protein occurs multiple
    times. 

    :arg1: TODO
    :returns: TODO

    """
    arg1.sort(key=lambda x: (float(x[2])))
    CreatedRange = range(0, len(arg1))
    for i in CreatedRange:
        CheckThisArea = range(int(arg1[i][2])-10, int(arg1[i][2]) + 45 )
        for item in arg1[i + 1 :]:
            if int(item[2]) in CheckThisArea:
                return 1
                break
            else:
                return 0
        break


def DuplicationSplitter(BlastListofList):
    """TODO: Docstring for DuplicationSplitter.

    :BlastListofList: TODO
    :returns: TODO

    """
    sorted(BlastListofList, key = lambda x: (x[2],x[4]))
    for item in BlastListofList:
        print item
    print '\n'

    def FindMultipleStarts(arg1):
        """TODO: Docstring for FindMultipleStarts.

        :arg1: TODO
        :returns: TODO

        """
        pass


def FileWriter(ListofBlastData):
    """TODO: Docstring for FileWriter.

    :ListofBlastData: TODO
    :returns: TODO

    """

    Scaffoldname = ListofBlastData[0][1].replace("|","_")
    ScaffoldFileName = str(Scaffoldname) + "blastTable.blast"


#    if os.path.exists(ScaffoldFileName):
#        os.remove(ScaffoldFileName)
    
    with open(ScaffoldFileName, 'a+') as f:
        for blast in ListofBlastData:
            Newitem = '\t'.join([str(item) for item in blast])
            f.write(Newitem)
            f.write('\n')

def SeqFileWriter(ListofBlastData, Sequence):
    """TODO: Docstring for FileWriter.

    :ListofBlastData: TODO
    :returns: TODO

    """
    print ListofBlastData
    Scaffoldname = ListofBlastData[0][1].replace("|","_")
    ScaffoldFileName = str(Scaffoldname) + "sequence.fasta"
    print Scaffoldname

    print Sequence
    #ISSUE RIGHT NOW IS WE ARE PRINTING MULTIPLE TIMES TO THE FILE. WE WANT ONE
    #SEQUENCE, and ONE HEADER
    
    with open(ScaffoldFileName, 'a+') as f:
        GeneName  = ListofBlastData[0][0]
        NewGeneName = ">" + str(GeneName)
        f.write(NewGeneName)
        f.write("\n")
        f.write(Sequence)
        f.write("\n")









def Usage():
    print "\n Application %s [options] -i <SeqFile> -b >blastfile> \n" \
        "-i     Input Genomic Scaffold File \n" \
        "-b     BlastFile with tab delinieated colulmn in specified format \n" \
        "-o     The output file you are going to write to. DO NOT INCLUDE.ending. All FILES will end with a ,out \n" \
        "'7 qseqid qstart qend sseqid sstart send length pident' | \n" \
        "sort -k 3,3n > BLASTINPUTFILE \n" % (sys.argv[0])

def Main():

    global oflag
    Iflag = None
    Bflag = None
    oflag = None
    Nucflag = None

    try:
        options, other_args = getopt.getopt(sys.argv[1:], "i:b:h:o:", ["help"])

    except getopt.GetoptError:
        print "There was a command parsing error"
        Usage()
        sys.exit(1)

    for option, value in options:
        if option == "-i":
            Iflag = value
        elif option == "-b":
            Bflag = value
        elif option == "-o":
            oflag = value
        elif option == "--help":
            Usage()
        else:
            print "Unhandeled options %s %s" % (options)

    if Iflag == None:
        print "Need a Genomic Scaffold input"
        Usage()
        exit(-1)
    elif Bflag == None:
        print "Need a blast file to referance"
        Usage()
        exit(-1)
    elif oflag == None:
        print "Need output"
        Usage()
        exit(-1)
    
    starttime = datetime.now()
    global Genome
    Genome = GenomeReader(Iflag)
    ReadBlast = BlastFilereader(Bflag)
    DictionaryOfScaff = Scafs2Dict(ReadBlast)
    GeneScafDict = BlastParser(ReadBlast, DictionaryOfScaff)
    CleanedDict = DictEditor(GeneScafDict)
    
    for Scaf, Gene in CleanedDict.iteritems():
        for GeneName, Blasthit in Gene.iteritems():
           GeneLenExonFiltered = PullingOutLowExonGenes(Blasthit)
           if GeneLenExonFiltered == 1:
               #NEEd to PUT IN DIRECTION SOMEWHERE?
               #FileWriter(Blasthit)
               GeneSeq = GeneObj(Blasthit)
               GeneSeq.FindGeneDirection()
               GeneSeq.RemoveOverlap()
               GeneSeq.SeqRetrieval(Genome)
               #SeqFileWriter(Blasthit, GeneSeq.Seqq)

           else:
               FindDuplicationStatus = FindDuplications(Blasthit)
               if FindDuplicationStatus == 0:
                   #FileWriter(Blasthit)
                   GeneSeq = GeneObj(Blasthit)
                   GeneSeq.FindGeneDirection()
                   GeneSeq.RemoveOverlap()
                   GeneSeq.SeqRetrieval(Genome)
                   #SeqFileWriter(Blasthit, GeneSeq.Seqq)

                    
                   #Call GENE OBJ
               elif FindDuplicationStatus == 1:
                   DuplicationSplitter(Blasthit)


    
    print datetime.now() - starttime






if __name__ == '__main__':
    Main()
