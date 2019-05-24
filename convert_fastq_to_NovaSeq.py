#!/usr/bin/env python2

__author__ = "Hannah Holland-Moritz"
__email__ = "hannah.hollandmoritz@colorado.edu"
__version__ = "0.0.1"

"""
Change MiSeq fastq files to match NovaSeq files
This program depends on biopython; to install use 

pip install biopython

learn more about biopython here: https://biopython.org/
"""

import argparse
from Bio import SeqIO
import gzip


def main():
    # Create the help messages and input requirements for the script

    parser = argparse.ArgumentParser(description=\
    '''
    Modifies SILVA taxonomy fasta file to be complient with the USEARCH sinTAX command.
    This program depends on biopython; to install it use
    
    pip install biopython
    ''',
    formatter_class=argparse.RawTextHelpFormatter)
    # what follows is some fancy code so /home/hannah/Documents/Fierer_lab/NovaSeqTestdada2that the required arguments are printed
    # above the optional arguments in the --help display.
    parser._action_groups.pop()
    req = parser.add_argument_group('required arguments')
    opt = parser.add_argument_group('optional arguments')

    # now the real arguments
    req.add_argument('-i', '--input_fp', required=True, type=str,
                     help=
     '''
     The input fastq file. 
     ''')

    req.add_argument('-o', '--output_fp',
        required=True, type=str,
                     help=
     '''
     The output file path.
     ''')

    opt.add_argument('-l', '--length', type=int,
                     help=
     '''
     The number of sequences in input file. This argument
     is not required but may speed up processing for long files
     if it is included. If you choose to use this argument,
     ***PLEASE MAKE SURE IT IS CORRECT***, otherwise, the program
     will take the time to read through your entire fastq file first
     before failing which will be a very sad day for you if your 
     fastq file is massive.
     ''')

    # Now the actual code

    # First extract the arguments  from the command line and save the input and
    # output file and the length, if it is given
    args = parser.parse_args()

    input_fp = args.input_fp
    output_fp = args.output_fp
    length = args.length

    # create a variable of correct length (if known) to hold modified sequences
    if length is not None:
        mod_records = [None] * length
        i = 0
    else:
        mod_records = []

    # Check if file is compressed or not
    iscompressed = check_file_compression(input_fp)
    
    if iscompressed:
        input_handle = gzip.open(input_fp, "rt")
    else:
        input_handle = input_fp


    print("Modifying sequence quality scores...")

    # read through the fastq file, for every record run the modify_fasta_header function
    for record in SeqIO.parse(input_handle, "fastq"):
        # extract header from record
        # print(record.letter_annotations["phred_quality"]) # debugging
        # print(record) # debugging

        # modify the quality information in NovaSeq format
        qual_out = modify_fastq_quality(list(record.seq),
                                        record.letter_annotations["phred_quality"])

        # overwrite the record quality score
        record.letter_annotations["phred_quality"] = qual_out

        # save the record to the modified records list;
        # if list is not pre-allocated in size, simply append the records
        # otherwise, incrementally replace value of "None" in the records list
        # with the record.
        if length is None:
            mod_records.append(record)
        else:
            mod_records[i] = record
            i += 1
        # print(mod_records) # debugging

    # If the input length doesn't match the actual number of sequences
    # exit with an error; Note, this is a stupid place to exit from the program,
    # since it exists after doing all the work, but if these two don't match,
    # the program will exit anyway, since writing the output file will fail.
    # This is just to make it fail with an error that makes sense.
    if length is not None and len(mod_records) != i:
        raise RuntimeError('It looks like the number of sequences you have specified '
                           'in the --length option does not match the actual number of '
                           'sequences in the input fasta file. Please recount the '
                           'sequences and try again.')


    print("Writing modified headers to file...")

    # Write modified records to file
    if iscompressed:
        with gzip.open(output_fp, "wb") as handle:
            SeqIO.write(mod_records, handle, "fastq")
    else:
        with open(output_fp, "w") as handle:
            SeqIO.write(mod_records, handle, "fastq")


# Helper functions #
def check_file_compression(input_fp):
    
    """
    check_file_compression(input_fp)
    This file uses the gzip "magic number" to check if a file is compressed or
    not and returns True if it is, or False if it is not.
    
    """
    
    gzip_magic_number = "1f8b"
    
    f = open(input_fp)
    
    iscompressed = f.read(2).encode("hex") == gzip_magic_number
    
    return iscompressed
    
    
def  modify_fastq_quality(fastq_seq, fastq_qual):

    """
    modify_fastq_quality(fasta_header)
    This function takes fastq sequence and fastq traditional numeric quality 
    score list as input and converts it to the new NovaSeq quality equivalents. 
    According to illumina 
    (https://www.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/novaseq-hiseq-q30-app-note-770-2017-010.pdf),
    quality scores <15 are assigned 12 ('-'), scores between 15 and 30 
    (inclusive) are assigned 23 ('8'), and scores >30 are assigned 37 ('F'). 
    Any null scores (i.e. base calls of N) are assigned 2 ('"'). 
    
    The function does the following steps:
        1. First checks for the presence of Ns in the sequence.
        2. If there are Ns, replaces the corresponding quality score with a 2, 
        otherwise continue to step 3
        3. Then translate scores <= 2 as 2, >2 and <15 as 12, >=15 and <31 as 
        23, and any score >30 (i.e. 31 and up) as 37.
    
    """
    # First checks for the presence of Ns in the sequence.
    if 'N' in fastq_seq:
    
       # If there are Ns, replaces the corresponding quality score with a 2
       for n, i in enumerate(fastq_seq):
           if i == 'N':
               fastq_qual[n] = 2
               
    # Translate scores           
    for n, i in enumerate(fastq_qual):
        if i <= 2:
            fastq_qual[n] = 2
        elif i < 15:
            fastq_qual[n] = 12
        elif i < 31:
            fastq_qual[n] = 23
        else:
            fastq_qual[n] = 37


    modified_fastq_qual = list(fastq_qual)

    return modified_fastq_qual


if __name__ == "__main__":
    main()
