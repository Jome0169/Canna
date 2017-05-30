import os
import sys
import getopt
from Bio.Seq import Seq
from datetime import datetime


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

def DictioanryOfScaffolds(ListofBlast):
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

#
#def PullOutIdenticalHits(DictOfBlastGenesandHits):
#    """Function to remove who hit scaffold in same location. Makes choice of
#    whcih to keep based off of e-score. Lower == Better
#
#    :DictOfBlastGenesandHits: TODO
#    :returns: TODO
#
#    """
#    #arg1.sort(key=lambda x: (float(x[2])))
#    ListofKeys = []
#
#    print DictOfBlastGenesandHits
#    for scaffold, gene in DictOfBlastGenesandHits.iteritems():
#        GeneNameSplit = gene.split("__")
#        if GeneNameSplit[1] not in ListofKeys:
#            ListofKeys.append(GeneNameSplit)
#        elif GeneNameSplit[1] in ListofKeys:
#            print "WE FOUND DUPLICATES PROBS"
#            print ListofKeys
#            print GeneNameSplit[1]
#
#
#
#
#        print scaffold
#        print gene
#        for gene2, locations in gene.iteritems():
#
#            print locations
#        print '\n' 
#            



def PullingOutLowExonGenes(arg1):
    """TODO: Docstring for PullingOutLowExonGenes.

    :arg1: Takes in list of blast hits that have less than 5 exons. Individuals
    with seq less than 800 Nts in length are thrown out. None type is returned
    in that case.  The idea here is that they might be psuedo genes or just have a few segments.
    However there are quite a few indivuduals with large hit boxes. So these
    will be pulled

    """


    if len(arg1) == 1 and (int(arg1[0][3]) - int(arg1[0][2])) <= 100:
        #Length Too Small
        return None
    elif len(arg1) == 1 and (int(arg1[0][3]) - int(arg1[0][2])) >= 600:
        #Probably Good Gene. return arg
        return arg1
    elif len(arg1) > 1: # If greater than length 1, begin more filtering 
        y = Filter2(arg1)
        if len(y) == 0:
            return None
        else:
            return y

def Filter2(arg2):
    """ If there is more than 1 exon we have to spend some time
    sorting the exons as well as looking for duplication. Takes in list,
    ensure total sum is greater than 800. Checks for overlap between the
    genes using protein query hit lovation (if both start with one)

    :arg2: List of list. Genes with less than 5 exons but more than 1. If
    total len is less than 800 Seq is deleted and Nonetype retruned
    :returns: None type or a list of list

    """
    NucLenSum = 0
    for item in arg2: # Calculate how much Nucleotide seq there is 
        SeqLen = (int(item[3]) - int(item[2]))
        NucLenSum += SeqLen 
    if NucLenSum >= 800: # if Greater than 800 nuceltodies (260 AA),proceed
        pass
    else:
        #Gene too short. Delete Argument. Will return None
        del[arg2[0:]]
    return arg2




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


def SplitDuplicateGenes(arg1):
    """TODO: Docstring for SplitDuplicateGenes.

    :arg1: TODO
    :returns: TODO

    """
    

    arg1.sort(key=lambda x: (float(x[2])))
    CopyList = arg1
    PossibleStart = []
    ProtStart = int(arg1[0][2])
    CheckThisArea = range(ProtStart - 5, ProtStart + 45) #Checks around prot strt
    PossibleStart.append(arg1[0])
    if len(arg1) % 2 == 0:
        print "EVEN SHOULD WORd"
        for item in arg1:
            print item
    else:
        print "ODD SHIT"
        for item in arg1:
            print item
    print '\n'
    #for item in arg1[1:]:
    #    if int(item[2]) in CheckThisArea:
    #        PossibleStart.append(item)
    #for item in PossibleStart:
    #    CopyList.remove(item)
    #if len(CopyList) == 0:
    #    pass
    #else:
    #    print arg1
       #Pairer(PossibleStart, CopyList)
       # print "STARTS"
       # print PossibleStart
       # print "POSS CONTINU"
       # print CopyList
       # print '\n'
   
    def Direction(arg1):
        """TODO: Docstring for Direction.

        :arg1: TODO
        :returns: TODO

        """
        Direction1 = ''
        if int(arg1[4]) < int(arg1[5]):
            Direction1 += '+'
        elif int(arg1[4]) > int(arg1[5]):
            Direction1 += '-'
        return Direction1
         

    def Pairer(list1, list2):
        """TODO: Docstring for Pairer.

        :list1: TODO
        :list2: TODO
        :returns: TODO

        """
        #if len(list1) == len(list2):
        for item in list1:
            Strand = Direction(item)
            if Strand == "+":
                Difference = []
                for thing in list2:
                    OtherStrand = Direction(thing)
                    if Strand == OtherStrand:
                        Diff = int(thing[5]) - int(item[4])
                        #print thing[5], item[4]
                        #print Diff
                        Difference.append(Diff)
                if Difference:
                    print Difference
                    K = min([n for n in Difference if n>0])
                    if K and int(K) <= 19000:
                        IndextoTake = Difference.index(K)
                        print item 
                        print list2[IndextoTake]
                        print '\n'
                else:
                    pass

#
#                    #TakeLowest Number
#                    #GetIndex 
#                    #Pring Lowest
#                print Difference
#
#
#            elif Strand== '-':
#                pass
#
#            
#


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
    ReadThatBlast = BlastFilereader(Bflag)
    DictionaryOfScaff = DictioanryOfScaffolds(ReadThatBlast)
    Z = BlastParser(ReadThatBlast, DictionaryOfScaff)
    for Scaffold, Gene in Z.iteritems():
        print Scaffold
        for Gene2, BlastLists in Gene.iteritems():
            print Gene2

        print '\n'
            #Bork = PullingOutLowExonGenes(BlastLists)
            #if Bork == None:
            #    pass
            #else:
            #    DuplicationStatus = FindDuplications(Bork)
            #    if DuplicationStatus == 0:
            #       CALL TO CLASS OBJ GENECLASS.py
            #        #pass
            #        #Z.SeqRetrieval(Genome)
            #    elif DuplicationStatus == 1:
            #        #SplitDuplicateGenes(Bork)
            #        pass


    print datetime.now() - starttime





if __name__ == '__main__':
    Main()


