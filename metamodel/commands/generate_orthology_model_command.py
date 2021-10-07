'''
Usage:
    metamodel generate_model [--help] [options] [--out=<PATH>]
                                       [--reaction_table=<PATH>] [--kegg_reference=<PATH>]

Required Arguments:
    out       	        Path to the directory you would like the output 
                        files to go (This will be created)
    reaction_table      A two-column table showing associations between 
                        the gene names of the proteome and a list of kegg
                        reaction IDs (delimited with a comma). This must 
                        be generate outside of this package, as we do not
                        currently have the capacity to generate orthology
                        predictions here
    kegg_reference      The location of the kegg database in a yaml format.
                        This can be generated using the generate_kegg_reference
                        command
                        
'metamodel generate_proteome' will generate a draft metabolic model that is compatible with
psamm based on eggnog predictions of the putative proteom identified from the generate_proteome
command.
'''

from docopt import docopt

if __name__ == '__main__':
    print(docopt(__doc__))
