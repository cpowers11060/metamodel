'''
Usage:
    metamodel generate_proteome [--help] [options] [--NCBI=<ACCESSION>] [--out=<PATH>]
                                [--assembly=<PATH>] [--pident=<INT>] [--pcover=<INT>]

Required Arguments:
    NCBI      A full NCBI accession. e.g. GCF_002020725.1_ASM202072v1
    out       Path to the directory you would like the output files to go (This will be created)
    assembly  Path to the fasta files containing all proteins from the assembly

Optional Arguments:
    pident    Specify a percent identity threshold for a protein to be considered
                a match to the reference [default: 50]
    pcover    Specify a percent cover threshold for a protein to be considered a
                match to the reference [default: 50]

'metamodel generate_proteome' downloads a specified NCBI taxa and uses a
fasta file from a metatranscriptome or a metagenome in order to generate
a putative proteome based on proteins present as good hits from the
reference proteome.
'''

from docopt import docopt

def main():
    print(docopt(__doc__))

if __name__ == '__main__':
    main()
