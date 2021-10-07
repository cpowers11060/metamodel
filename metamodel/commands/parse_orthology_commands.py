'''
Usage:
    metamodel parse_orthology [--help] [options] [--orthology=<PATH>] [--out=<PATH>]
                                [--strategy=<STR>]

Required Arguments:
    orthology   The location of the orthology table
    out         Path to the directory you would like the output files to go (This will be created)
    strategy    Strategy for the generation of the orthology table, such as Eggnog. Default eggnog

'metamodel parse_orthology' takes output from an orthology table and generates a new
table with gene to reaction relationships.
'''

from docopt import docopt

if __name__ == '__main__':
    print(docopt(__doc__))
