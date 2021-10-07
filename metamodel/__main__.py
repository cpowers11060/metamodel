#! /usr/bin/env python
'''
Usage:
    metamodel [--version] [--help] <command> [<args>...]

Options:
    -h, --help          Generate Help Screen
    -v, --version       Get Version Number

General Commands:
    generate_proteome   Generates a draft proteome of a putative
                        organism based on pairwise alignments 
                        between a reference file and the input
                        assembly. Outputs the following:
                        - fasta file of proteome
                        - configuration files for Circos
                        - run summary information
    generate_model      Generates a draft metabolic model based
                        on a specified draft model from the psamm
                        model collection. (Dependent on running
                        generate_proteome first). Generates the
                        following:
                        - A reciprocal best hits file between the 
                        proteome and the reference model proteome
                        - A new model that shows the functions 
                        represented in the reciprocal best hits
                        file
                        - A visualization of the draft metabolic
                        network
    parse_orthology     Takes an input orthology mapping and generates
                        a reaction:gene association table for the
                        generate_model script. Currently accepts:
                        - Eggnog tabular output

See 'metamodel <command> --help' for more information on a command
'''

import os
from docopt import docopt

import metamodel.run_psammotate
import metamodel.generate_model
import metamodel.generate_proteome
import metamodel.visualize_draft_model
import metamodel.parse_orthology

if __name__ == '__main__':
    args = docopt(__doc__,
                  version='',
                  options_first=True)
    argv = [args['<command>']] + args['<args>']

    if args['<command>'] == 'generate_proteome':
        import commands.generate_proteome_command
        args = docopt(commands.generate_proteome_command.__doc__,
                      argv=argv)
        if args['--out']:
            if os.path.exists(args['--out']):
                raise Exception('The path already exists! We do not want to overwrite data...')
        else:
            args['--out']=="./proteome_out"
        if args['--assembly']:
            if not os.path.exists(args['--assembly']):
                raise Exception('There is no file at this location! Assembly does not exist!')        
        else:
            raise Exception('You didnt specify an assembly, you should do that')
        if args['--NCBI'] is None:
            raise Exception('No accession was specified with the --NCBI flag. Exiting...')
        if args['--pident'] is None:
            args['--pident']=50
        if args['--pcover'] is None:
            args['--pcover']=50
        metamodel.generate_proteome.create_prot(args['--NCBI'],
                                                 args['--out'],
                                                 args['--assembly'],
                                                 args['--pident'],
                                                 args['--pcover'])

    elif args['<command>'] == 'generate_model':
        import commands.generate_orthology_model_command
        args = docopt(commands.generate_orthology_model_command.__doc__,
                      argv=argv)
        if args['--out']:
            if os.path.exists(args['--out']):
                raise Exception('The path already exists! We do not want to overwrite data...')
        else:
            args['--out']=="./model_out"
        if args['--reaction_table']:
            if not os.path.exists(args['--reaction_table']):
                raise Exception('There is no file at this location! Reaction Table does not exist!')
        else:
            raise Exception('You didnt specify a Reaction Table, you should do that')
        if args['--kegg_reference']:
            metamodel.generate_model.create_model_kegg(args['--out'],args['--reaction_table'],
                                              args['--kegg_reference'])
        else:
            metamodel.generate_model.create_model_api(args['--out'],args['--reaction_table'])
    elif args['<command>'] == 'parse_orthology':
        import commands.parse_orthology_commands
        args = docopt(commands.parse_orthology_commands.__doc__,
                      argv=argv)
        if args['--out']:
            if os.path.exists(args['--out']):
                raise Exception('The path already exists! We do not want to overwrite data...')
        else:
            args['--out']=="./ortho_out"
        if args['--orthology']:
            if not os.path.exists(args['--orthology']):
                raise Exception('Why would you not specify an orthology table. Do you not want this to work?')
        if not args['--strategy']:
            args['--strategy']='eggnog'
        if args['--strategy']=='eggnog':
            metamodel.parse_orthology.parse_eggnog(args['--out'],args['--orthology'])
        else:
            raise Exception('Please specify a valid strategy. Try --help if you don\'t know what I\'m talking about')
