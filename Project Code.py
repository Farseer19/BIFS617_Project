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

# function for reverse complement
def reverse_complement(seq):
    # replace nucleotides with appropriate complement
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    # reverse the sequence
    return seq.translate(complement)[::-1]


# Finds ORFs in a given frame 
def find_orfs_in_frame(seq, min_len=50, reverse = False):
    # List to store found ORFs
    orfs = []

    # Loop through all three frame offsets (0, 1, or 2)
    for offset in range(3):
        # Start reading at the given frame offset 
        i = offset

        # Loop through the sequence in steps of 3
        while i <= len(seq) - 3:
            # Extract the codon at current position
            codon = seq[i:i+3]

            # If it's a start codon, look for the stop codon
            if codon == 'ATG':
                # Mark the start of the ORF
                start = i 
                
                # Use j to identify stop codon while keeping start = i
                j = i      

                # Find a stop codon or reach end of sequence
                while j <= len(seq) - 3:
                    next_codon = seq[j:j+3]

                    # If a stop codon is found, check ORF length
                    if next_codon in ['TAA', 'TAG', 'TGA']:
                        
                        # ORF ends after the stop codon
                        end = j + 3  
                        
                        # add in location of start and stop
                        orf_seq = seq[start:end]

                        # Only record ORFs that meet the minimum length
                        if len(orf_seq) >= min_len:
                            
                            # Frame number depending on which strand
                            frame = offset + 4 if reverse else offset + 1
                            
                            # Position number depending on which strand
                            pos = -(len(seq) - start) if reverse else start + 1
                            
                            # Add in orf if long enough  
                            orfs.append((frame, pos, orf_seq))  
                        
                        break  

                    # If orf is too short, keep looking at next codon
                    j += 3

                # Continue search with next start codon
                i = j + 3
            else:
                # Not a start codon, move to next codon
                i += 3  

    # Return list of ORFs found in this frame
    return orfs
