f"Algorithm to check each reading frame for start and stop codons, storing their positions in the original string"


f"Parse input FASTA files that may have errors  "


def input_FASTA ():
    f"Dictionary that stores header information and corresponding full sequences in all caps"
    return_FASTA = dict()
    with open(r"C:\Users\halld\Documents\input.txt") as f:
        f"Sets current FASTA entry where sequence data is being added to"
        curr_FASTA = str()

        for line in f:
            line=line.strip()
            if line:
                if line[0] == '>':
                    curr_FASTA = line
                    return_FASTA[curr_FASTA] = str()
                else:
                    return_FASTA[curr_FASTA] += line.upper()
    return return_FASTA



print(input_FASTA())