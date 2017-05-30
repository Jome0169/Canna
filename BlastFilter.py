import os
import sys
import getopt
from Bio.Seq import Seq
from datetime import datetime




def BlastFilereader(BLastfile):
    # BlastFilereader will return a list of list from a blast file with the
    #above commented 7 header feature from a BLAST table."""
    HitsNStuff = []
    with open(BLastfile, "r") as f:
        for line in f:
           HitsNStuff.append(line.strip("\n").split('\t'))
        
    HitsNStuff.sort(key=lambda x: x[1])
    return HitsNStuff


def Scaffoldreturn(arg1):
    """TODO: Docstring for BlastEditer.

    :arg1: TODO
    :returns: TODO

    """
    ListForScaffolds = []
    for item in arg1:
        if item[1] not in ListForScaffolds:
            ListForScaffolds.append(item[1])
    return ListForScaffolds


def WriteFiles(ScaffoldNames, Blastfile):
    """TODO: Docstring for WriteFiles.

    :ScaffoldNames: TODO
    :ListofScaffolds: TODO
    :returns: TODO

    """

    for item in ScaffoldNames:
        NewFileName = str(item).replace("|", '_') + "_scaffold.blast"
        LinesToadd = []

        with open(Blastfile, "r") as f:
            for line in f:
                Newline = (line.strip("\n").split('\t'))
                if item in Newline[1]:
                    LinesToadd.append(line)
        with open(NewFileName, 'a+') as Z:
            for item in LinesToadd:
                Z.write(item)



def Usage():
    print "\n Application %s [options] -i <SeqFile> -b >blastfile> \n" \
        "-i     Input Genomic Scaffold File \n" \
        "-b     BlastFile with tab delinieated colulmn in specified format \n" \
        "-o     The output file you are going to write to. DO NOT INCLUDE.ending. All FILES will end with a ,out \n" \
        "'7 qseqid qstart qend sseqid sstart send length pident' | \n" \
        "sort -k 3,3n > BLASTINPUTFILE \n" % (sys.argv[0])


def Main():

    Bflag = None

    try:
        options, other_args = getopt.getopt(sys.argv[1:], "b:h:", ["help"])

    except getopt.GetoptError:
        print "There was a command parsing error"
        Usage()
        sys.exit(1)

    for option, value in options:
        if option == "-b":
            Bflag = value
        elif option == "--help":
            Usage()
        else:
            print "Unhandeled options %s %s" % (options)

    if Bflag == None:
        print "Need a blast file to referance"
        Usage()
        exit(-1)

    Raw_blastordered = BlastFilereader(Bflag)
    Scaffolds = Scaffoldreturn(Raw_blastordered)
    WriteFiles(Scaffolds, Bflag)

 




if __name__ == '__main__':
    Main()
