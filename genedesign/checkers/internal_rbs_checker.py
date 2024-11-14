from genedesign.seq_utils.reverse_complement import reverse_complement


class InternalRBSChecker:

    def __init__(self):
        self.shines = None

    def initiate(self):
        self.shines = ["AGGAGG", "GGAGG", "AGGAG"]

    def mismatches(self, shine, seq):
        if len(shine) > len(seq):
            return len(shine)
        minDiff = len(shine)
        for i in range(0, 1 + len(seq) - len(shine)):
            count = 0
            for j in range(len(shine)):
                #print(seq[i+j] + shine[j])
                if not seq[i+j] == shine[j]:
                    count += 1
            minDiff = min(count, minDiff)
        return minDiff
    
    def run(self, sequence):

        rc = reverse_complement(sequence)

        seq = sequence + "           " + rc

        for i in range(0, len(seq) - 3):
            if seq[i: i + 3] == "ATG":
                for shine in self.shines:
                    #print(seq[max(0, i - len(shine) - 8):i - 4])
                    if self.mismatches(shine, seq[max(0, i - len(shine) - 8):i - 5] ) <= 1:
                        return False,  seq[max(0, i - len(shine) - 8):i + 3]
        
        return True, None

#print(r.mismatches("AGGAGG", "AGGATT"))
#print(r.run("AGGAGGAAAATG"))