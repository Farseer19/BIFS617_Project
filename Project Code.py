f"Algorithm to check each reading frame for start and stop codons, storing their positions in the original string"


f"Parse input FASTA files that may have errors  "

#Author: Duncan Hall
#Function for handling input of FASTA file
def input_FASTA ():
    f"Dictionary that stores header information and corresponding full sequences in all caps"
    return_FASTA = dict()
    #Opening stored FASTA file
    with open(r"input.txt") as f:
        f"Sets current FASTA entry where sequence data is being added to"
        curr_FASTA = str()
        #Parse through each line in file
        for line in f:
            #Strip each line in file
            line=line.strip()
            #Checks if there is anything in this "line"
            if line:
                #If the line is the header, the key in the dictionary will be set that with a string as the value to be extended as sequence lines are read
                if line[0] == '>':
                    curr_FASTA = line
                    return_FASTA[curr_FASTA] = str()
                #Sequence is extended while there is a line and that line is not a header
                else:
                    return_FASTA[curr_FASTA] += line.upper()
    #Return the FASTA dictionary
    return return_FASTA




# Author: Matt Kubit
# function for reverse complement
def reverse_complement(seq):
    # replace nucleotides with appropriate complement
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    # reverse the sequence
    return seq.translate(complement)[::-1]

# Author: Matt Kubit
# Finds ORFs in a given frame 
def find_orfs_in_frame(seq, min_len, reverse):
    # List to store found ORFs
    orfs = []

    if reverse:
        seq = reverse_complement(seq)
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
                            pos = -(start+1) if reverse else start + 1
                            
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

#Author: Duncan Hall
#Function to call input_FASTA function and find_orfs_in_frame in conjunction to create the output. Only variable is min_len so the desired minimum length can be passed to the find_orfs_in_frame call.
def orfs_from_input(min_len):
    #Call input_FASTA to retrieve sequence data
    input_dict = input_FASTA()
    #Create output_dict to store output data
    output_dict = {}
    f"Look through each entry from the input_FASTA function and process to find ORFs"
    for entry in input_dict:
        #Calls find_orfs_in_frame in the forward direction and stores the data in output_dict
        output_dict[entry] = find_orfs_in_frame(input_dict[entry],min_len,True)
        #Calls find_orfs_in_frame in the reverse direction and stores the data in output_dict
        output_dict[entry] += find_orfs_in_frame(input_dict[entry],min_len,False)
    
    print(output_dict)
    f"Print out in viewable format"
    for entry in output_dict:
        for orf in output_dict[entry]:
            print(entry, "|FRAME=",orf[0],"POS=",orf[1],"LEN=",len(orf[2]))
            output_seq = " ".join([orf[2][i:i+3] for i in range(0, len(orf[2]), 3)])
            print(output_seq)

f"Printing with no minimum"
orfs_from_input(0)
