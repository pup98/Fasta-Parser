import wget                     
from Bio import SeqIO
import sys

# Download the multifasta file

def download_file(url):
    try: 
        file = wget.download(url)
        return file
    except:
        return('URL is broken!')

# Let's save at all the fasta entries it has in a variable

def read_fasta(file):
    with open(file) as my_file:
        output_fasta = my_file.read()
     
    '''Now we will convert the output fasta string containing all the sequences to a list, 
    to help accessing the sequence that we want to process'''

    sequence_list = output_fasta.split(">")
    return sequence_list

''' We will pick second element of sequence_list to perform mutation and realted processing and
append the sequence under study to a new file. Additionally, remove the working sequence from multifasta file
and process further'''

def extract_sequence(sequence_list):

    '''Since we splitted the output_fasta by '>', this charcter will be absent for the sequences in the list!
    . Therefore, below line to convert it back to a fasta format''' 
    extracted_seq = '>'+ sequence_list[1]

    # Save the working sequence in a new file
    with open('extracted_seq.fasta','w') as outfile:
        outfile.write(extracted_seq)
    
    # Let's look at specifics of the sequence under study and return: id, sequence
    with open('extracted_seq.fasta') as file:
        for seq_record in SeqIO.parse(file, "fasta"):
            
            removal_record_id = seq_record.id
            removal_seq = seq_record.seq
    
    '''Parse the parent file and append the rest seq.id's except the working seq to a new file
     This will append the new file with out our sequence of concern'''

    with open('mutated_multifasta.fn','w') as multifasta_file, open('extracted_seq.fasta','r') as output_file:
        remove = output_file.read() 
        
        for seq in SeqIO.parse('Vieuvirus.fn', 'fasta'):
            if seq.id not in remove:
                SeqIO.write(seq, multifasta_file, "fasta")
            
    '''removal_seq is a Bio.seq object, let's convert it into string for ease of processing'''
    removal_seq = str(removal_seq)

    return removal_seq, removal_record_id


# Function for reverse complement

def find_reverse_complement(dna_sequence):
    try: 
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        reverse_complement_sequence = ''.join(complement_map[base] for base in dna_sequence[::-1])
        return reverse_complement_sequence
    except:
        return('DNA sequence has garbage characters! check the sequence again')


# Mutating the reverse complement sequence at random position

# First lets convert the str object to list, for ease of processing
def mutate_addtomultifasta_seq(reverse_complement_sequence, removal_record_id):
    reverse_complement_sequence = list(reverse_complement_sequence)
    reverse_complement_sequence[2] = "C"
    reverse_mutated = ''.join(base for base in reverse_complement_sequence)
    
    # Convert the reverse complement to a Seq object, before inserting to multifasta file 
    from Bio.Seq import Seq
    mutated_seq = Seq(reverse_mutated)

    # Convert the sequence into seq record and append to a new file
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(mutated_seq, removal_record_id)
    SeqIO.write(record, "mutated_seq.fasta", "fasta")

    '''Updating the multifasta 'mutated_multifasta.fn' file by adding 'mutated_seq.fasta'  ....'''

    with open('mutated_multifasta.fn', "a") as f, open('mutated_seq.fasta', "r") as r:
        for seq in SeqIO.parse(r, 'fasta'): 
                SeqIO.write(seq, f, "fasta")


# Design and Update the bed file
# Annotate the bed file accordingly as to what modifications we did!

def create_bed_file(file, record_id):
    with open(file) as multifasta_file, open('updated_sequences.bed','w') as output_bed_file:
        for seq in SeqIO.parse(multifasta_file, 'fasta'):
            if seq.id == record_id:
                output_bed_file.write('{}\t0\t{}\t{}\n'.format(seq.id, len(seq), 'Mutated sequence'))
            else: 
                output_bed_file.write('{}\t0\t{}\t\n'.format(seq.id, len(seq)))


if __name__ == "__main__":
    url = 'https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/Vieuvirus.fn'
    file = download_file(url)

    sequence_list = read_fasta(file)

    removal_seq, removal_record_id = extract_sequence(sequence_list)

    reverse_complement_sequence = find_reverse_complement(removal_seq)

    mutate_addtomultifasta_seq(reverse_complement_sequence, removal_record_id)
    create_bed_file(file, removal_record_id)