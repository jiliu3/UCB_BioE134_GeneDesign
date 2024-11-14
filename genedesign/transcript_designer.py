from genedesign.seq_utils.hairpin_counter import hairpin_counter

from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript

from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker

from typing import List
import random
import re
from decimal import Decimal


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using the highest CAI codon for each amino acid.
    """

    codonDict = {}

    def initialize_codonDict(self):
        if bool(self.codonDict):
           return 
        for line in open("genedesign\data\codon_usage.txt"):
            fields = re.split(" |\t", line)
            aa = fields[1]
            freq = int(Decimal(fields[2]) * 100)
            codon = fields[0]
            if self.codonDict.get(aa) == None:
                self.codonDict[aa] = [codon] * freq
            else:
                self.codonDict[aa] += [codon] * freq  


    def __init__(self):
   
        self.aminoAcidToCodon = {}
        self.rbsChooser = None
        self.forbiddenSequenceChecker = None
        self.promoterChecker = None
        self.rbsChecker = None
        self.forbidden = []
        self.codonDict = {}
        self.cc = None
    
            
    def initiate(self) -> None:
        """

        Initializes the codon table and the RBS chooser.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()

        self.forbiddenSequenceChecker = ForbiddenSequenceChecker()
        self.forbiddenSequenceChecker.initiate()

        self.rbsChecker = InternalRBSChecker()
        self.rbsChecker.initiate()

        self.cc = CodonChecker()
        self.cc.initiate()
        

        random.seed(2)

        self.initialize_codonDict()

        # Codon table with highest CAI codon for each amino acid (for E. coli)
        self.aminoAcidToCodon = {
            'A': "GCG", 'C': "TGC", 'D': "GAT", 'E': "GAA", 'F': "TTC",
            'G': "GGT", 'H': "CAC", 'I': "ATC", 'K': "AAA", 'L': "CTG",
            'M': "ATG", 'N': "AAC", 'P': "CCG", 'Q': "CAG", 'R': "CGT",
            'S': "TCT", 'T': "ACC", 'V': "GTT", 'W': "TGG", 'Y': "TAC"
        }

    def window(self, peptide: str, start: int, end:int, preamble: List[str]):
        repeats = 200
        choices = {}

        for i in range(repeats):
            temp = []
            for aa in peptide[start:end]:
                size = len(self.codonDict[aa])
                temp.append(self.codonDict[aa][random.randint(0, size - 1)])
            
            whole = "".join(temp) + "".join(preamble)
            #print(temp)
            if self.forbiddenSequenceChecker.run(whole) and self.promoterChecker.run(whole) and self.rbsChecker.run(whole):
                hCount = hairpin_counter(whole)[0]
                if hCount == 0:
                    result = self.cc.run(temp)
                    choices[result[3] + 0.5 * result[2]] = temp
        
        
        if bool(choices):
            #print(max(choices.keys()))
            return choices[max(choices.keys())]
        else:
            return temp
        
                    
        
        





    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA and selects an RBS.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """

        #create starting six amino acids
        codons = ""
        temp = ""
        codons = []

        peptide = peptide + "*"
        
        codons += self.window(peptide, 0, 3, "")
        
        index = 3

        end = len(peptide) - len(peptide)%3 - 6
        for i in range(3, end, 3):
            preamble = codons[index - 3:]
            temp = self.window(peptide, index, index + 9, preamble)
            codons += temp[:3]
            index += 3

        preamble = codons[index - 3:]
        codons += self.window(peptide, end, len(peptide), preamble)


        # Build the CDS from the codons
        cds = ''.join(codons)

        # Choose an RBS
        selectedRBS = self.rbsChooser.run(cds, ignores)

        # Return the Transcript object
        return Transcript(selectedRBS, peptide, codons)

if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    cc = CodonChecker()
    cc.initiate()

    # Print out the transcript information
    print(transcript)
    

