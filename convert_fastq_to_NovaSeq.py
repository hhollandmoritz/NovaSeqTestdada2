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
import re


def main():
    # Create the help messages and input requirements for the script

    parser = argparse.ArgumentParser(description=\
    '''
    Modifies SILVA taxonomy fasta file to be complient with the USEARCH sinTAX command.
    This program depends on biopython; to install it use
    
    pip install biopython
    ''',
    formatter_class=argparse.RawTextHelpFormatter)
    # what follows is some fancy code so that the required arguments are printed
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

    # Get the taxonomy dictionary (before the loop so that it only has to be created once)
    tax_dict = create_taxonomy_dictionary(map_fp)
    print("Taxa dictionary created.")

    print("Modifying sequence headers...")

    # read through the fastq file, for every record run the modify_fasta_header function
    for record in SeqIO.parse(input_fp, "fasta"):
        # extract header from record
        header=record.description

        # modify the header to match the sinTAX format, add ";" to sequence id
        seq_id_out, header_out = modify_fasta_header(header, tax_dict)

        # overwrite the record description and id with the new header and new sequence id
        record.description = header_out
        record.id = seq_id_out

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
    with open(output_fp, "w") as handle:
        SeqIO.write(mod_records, handle, "fasta")
        #record.description = header_out


# Helper functions #


def create_quality_dictionary():
    """
    create_quality_dictionary()
    this function creates a dictionary of quality scores where the keys
    are traditional illumna fastq quality symbols and the values are the 
    new NovaSeq quality equivalents. According to illumina 
    (https://www.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/novaseq-hiseq-q30-app-note-770-2017-010.pdf),
    quality scores <15 are assigned 12 ('-'), scores between 15 and 30 
    (inclusive) are assigned 23 ('8'), and scores >30 are assigned 37 ('F'). 
    Any null scores (i.e. base calls of N) are assigned 2 ('"'). As our pipeline
    currently filters out Ns before processing, I do not bother assining 2 to
    no-calls.
    """

    # initialize dictionary
    qual_dict = {"!" : "-",
                "\"" : "-",
                "#" : "-",
                "$" : "-",
                "%" : "-",
                "&" : "-",
                "'" : "-",
                "(" : "-",
                ")" : "-",
                "*" : "-",
                "+" : "-",
                "," : "-",
                "-" : "-",
                "." : "-",
                "/" : "-",
                "0" : "8",
                "1" : "8",
                "2" : "8",
                "3" : "8",
                "4" : "8",
                "5" : "8",
                "6" : "8",
                "7" : "8",
                "8" : "8",
                "9" : "8",
                ":" : "8",
                ";" : "8",
                "<" : "8",
                "=" : "8",
                ">" : "8",
                "?" : "8",
                "@" : "F",
                "A" : "F",
                "B" : "F",
                "C" : "F",
                "D" : "F",
                "E" : "F",
                "F" : "F",
                "G" : "F",
                "H" : "F",
                "I" : "F"
                }

    return qual_dict


def multiple_replace(text, adict):

    """
    multiple_replace(text, adict)
    takes a text string and a dictionary of replacements
    allows multiple replacement of characters so that those characters
    only need to be specified in a dictionary once
    depends on the re module
    for more information on this function, see
    https://www.safaribooksonline.com/library/view/python-cookbook-2nd/0596007973/ch01s19.html
    """
    rx = re.compile('|'.join(map(re.escape, adict)))

    def one_xlat(match):
        return adict[match.group(0)]

    final_replace = rx.sub(one_xlat, text)

    return final_replace


def match_taxa(tax_dict, split_taxonomy):

    """
    match_taxa(tax_dict, split_taxonomy)
    this function will take the taxonomy dictionary and the
    split taxa string as input and assign it the proper
    d:/p:/c:/o:/f:/g:. If a split is the last one in the
    list, it will be given the species name. Non-conforming
    taxonomy levels are removed from the list. It will return a
    list of the taxonomy matches modified to include  d:/p:/c:/o:/f:/g:/s:
    #
    Note:
    this step would probably be better if vectorized, but
    I don't know how to do that easily, plus I want to keep
    the script dependencies small. So it is in a loop instead.
    """

    for taxa_level in list(range(len(split_taxonomy))):

        # save the current taxonomy level as a new variable
        current_tax_level = split_taxonomy[taxa_level]

        # then if the level is above species, look up it's abbreviation
        if taxa_level is not (len(split_taxonomy) - 1):
            # for abbreviations that conform to the sinTax style,
            # overwrite the split_taxonomy list with the new
            # style, if they are non-conforming, mark them by
            # labeling them "non-conforming".

            if tax_dict[current_tax_level] not in ["non-conforming:"]:
                split_taxonomy[taxa_level] = tax_dict[current_tax_level] + current_tax_level
            else:
                split_taxonomy[taxa_level] = "non-conforming"

        # if a taxa is the species level (i.e. the last level in the string)
        #  put an "s:" before it.
        else:
            split_taxonomy[taxa_level] = "s:" + current_tax_level

    # remove non-conforming taxonomy levels from the list
    if "non-conforming" in split_taxonomy:
        while "non-conforming" in split_taxonomy:
            split_taxonomy.remove("non-conforming")
        modified_split_taxonomy = split_taxonomy
    else:
        modified_split_taxonomy = split_taxonomy

    return modified_split_taxonomy


def modify_fasta_header(fasta_header, tax_dict):

    """
    modify_fasta_headers(fasta_header)
    This function takes a fasta header string as input and a taxonomy dictionary and
    1) splits the string into the sequence id and the taxonomy string
    2) splits the string based on semicolon placement
    3) places d:, p:, c:, o:, f:, g:, and s: at the appropriate places before each
    taxonomy level with the use of a helper function, match_taxa()
    4) removes special characters from the taxonomy strings by replacing them with underscores
    5) rebuilds the string and appends it to the sequence ids
    6) returns the sequence id with a ";" after it and the modified taxonomy header
    """

    # First split the header string on the first space and save the
    # sequence ID and taxonomy [second half]
    taxonomy_string = fasta_header.split(" ", 1)[1]
    sequence_id = fasta_header.split(" ", 1)[0]

    # Then split the string based on semicolons
    split_taxonomy = taxonomy_string.split(";")

    # for each split_taxonomy, check which taxa level it is, and rename it with the proper
    # d:/p:/c:/o:/f:/g: heading from the taxa_dict (species are not included in that list)
    # so the last split in the list that does not match any taxa level
    # will be given the species name.
    modified_taxa_list = match_taxa(tax_dict, split_taxonomy)

    # Then replace special characters with underscoes
    replacement_dict = {
        " " : "_",
        "(" : "",
        ")" : "",
        "/" : "",
        "," : "_"
    }

    modified_taxa_list_no_spc_char = []

    for taxa in range(len(modified_taxa_list)):
        #print(taxa)
        taxa_name = modified_taxa_list[taxa]
        #print(taxa_name)
        taxa_str_no_spc_char = multiple_replace(taxa_name, replacement_dict)
        #print(taxa_str_no_spc_char)
        modified_taxa_list[taxa] = taxa_str_no_spc_char
        #print(modified_taxa_list[taxa])

    # Then rebuild the string from the list by joining taxa list
    # together with a comma as the separator; also add "tax=" to the beginning
    full_taxonomy_string = "tax=" + ",".join(modified_taxa_list)

    # now append the modified string to the sequence id

    modified_header = full_taxonomy_string + ";"
    modified_seq_id = sequence_id + ";"

    return modified_seq_id, modified_header


if __name__ == "__main__":
    main()
