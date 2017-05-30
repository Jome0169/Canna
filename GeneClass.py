import os
import sys
import getopt
from Bio.Seq import Seq
from datetime import datetime


class GeneObj(object):

    """Docstring for Gene. """

    def __init__(self, BlastLoc):
        """TODO: to be defined1.

        :BlastHits: TODO

        """
        self.BlastLoc = BlastLoc
        self.Orientatin = ''
        self.Seqq = ''


    def RevComp(self):
        """
        Uses biopython to take reverse complement of nucleotoide sequence if it
        was found to be negative in previous steps
        """
        Y = Seq(self.Seqq).reverse_complement()
        self.Seqq = str(Y)


    def FindGeneDirection(self):
        """TODO: Docstring for FindGeneDirection.

        :arg1: TODO
        :returns: TODO

        """

        arg1 = self.BlastLoc
        Negative = 0
        Posotive = 0
        CheckForNestedList = any(isinstance(i, list) for i in arg1) #checks if list nested
        
        if CheckForNestedList == True:
            for item in arg1:
                if int(item[4]) < int(item[5]):
                    Posotive += 1 
                elif int(item[4]) > int(item[5]): #space 4 laerger == NEgative
                    Negative += 1
            if Posotive != 0 and Negative == 0:
                self.Orientatin = '+'
            elif Negative != 0 and Posotive == 0:
                self.Orientatin = '-'
            else:
                return None
        
        elif CheckForNestedList == False:
            if int(arg1[4]) < int(arg1[5]):
                    Posotive += 1 
            elif int(arg1[4]) > int(arg1[5]): #space 4 laerger == NEgative
                    Negative += 1
            if Posotive != 0 and Negative == 0:
                self.Orientatin = '+'
            elif Negative != 0 and Posotive == 0:
                self.Orientatin = '-'
            else:
                return None

    def RemoveOverlap(self):
        """TODO: Docstring for RemoveOverlap.
        :returns: TODO

        """

        def ListEdit(Difference, List, Orien):
            """TODO: Finds the overlap between these list and edits the query
            location and BLAST hit location based off of overlap. 

            [['c27027_g2_i3_1__gi|703127562|ref|XM_010105561.1|__5e-15+',
            '001035F|quiver', '13', '860', '6673', '5833', '90.36', '0.0'],
            ['c27027_g2_i3_1__gi|703127562|ref|XM_010105561.1|__5e-15+',
            '001035F|quiver', '857', '888', '5729', '5698', '100.00', '3e-07']]

            for instance there is query overlap above of 860 and 857. We need
            to change this so 857 Instead picks up at the appropriate location.

            :Difference: TODO
            :List: TODO
            :returns: TODO

            """
            if Orien == '+':
                ModdedList = List
                ModdedList[2] = int(ModdedList[2]) + Difference
                ModdedList[4] = int(ModdedList[4]) + Difference
                return ModdedList
            elif Orien == '-':
                ModdedList = List
                ModdedList[2] = int(ModdedList[2]) + Difference
                ModdedList[4] = int(ModdedList[4]) - Difference
                return ModdedList


        arg1 = self.BlastLoc
        Direction = self.Orientatin
        RangeToCheck = xrange(0, len(arg1))
        for i in RangeToCheck:
            EndOfHit = int(arg1[i][3])
            for item in arg1[i + 1:]:
                if EndOfHit >= int(item[2]):
                    
                    FoundDifference = (EndOfHit - int(item[2])) + 1
                    ListEdit(FoundDifference, item, Direction)
                    break
                else:
                    pass
        self.BlastLoc = arg1
           
        

    def SeqRetrieval(self, GenomeSeq):
        """TODO: Docstring for SeqRetrieval.

        :arg1: TODO
        :returns: TODO

        """
        Direction = self.Orientatin
        BlastHits = self.BlastLoc
        SeqTaken = self.Seqq
        if Direction == '+':
            for item in BlastHits:
                Seq2 = GenomeSeq[item[1]][int(item[4]) -1 :int(item[5]) + 1 ]
                SeqTaken += Seq2
            #print SeqTaken
            #print Direction
            self.Seqq = SeqTaken
        elif Direction == '-':
            BlastHits.sort(key=lambda x: (float(x[4])))
            for item in BlastHits:
                Seq2 = GenomeSeq[item[1]][int(item[5]) -1 :int(item[4]) + 1 ]
                SeqTaken += Seq2
            self.RevComp()
            #print Direction
            #print SeqTaken

            self.Seqq = SeqTaken


