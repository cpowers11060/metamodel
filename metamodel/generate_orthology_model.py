import re
import os
from psamm.datasource.native import NativeModel, ModelReader, ModelWriter
from psamm.expression import boolean
from collections import defaultdict

def create_model(out, rxn, kegg):
    # Initialize some stuff that store the reaction data
    # Stores the Reaction:Gene associations
    rxn_mapping=defaultdict(list)
    with open(rxn, "r") as infile:
        for line in infile:
            line=line.rstrip()
            listall=re.split("\t", line)
            genes=re.split(":", listall[1])
            rxn_mapping[listall[0]]=genes

    # First read in the entire kegg model and define the empty output
    # native model
    mr = ModelReader.reader_from_path('{}'.format(kegg))
    nm = mr.create_model()
    out_nm = NativeModel()
    mm = nm.create_metabolic_model()

    # Initialize some stuff to store compounds, ec, reactions, etc.
    compound_list = []
    in_kegg_ec = set()

    # Loop through reactions in the kegg model
    for reaction in nm.reactions:
        for ec in reaction.properties['enzymes']:
            in_kegg_ec.add(ec)
        if reaction.id in rxn_mapping.keys():
            out_nm.reactions.add_entry(reaction)
            for compound in reaction.equation.compounds:
                compound_list.append(compound[0].name)

    # Create a model writer object that will write the new native
    # model out
    out_mw = ModelWriter()

    with open('./{}/gene-association.tsv'.format(out), mode='w') as outfile:
        outfile.write('id\tgenes\n')
        for reaction in out_nm.reactions:
            if len(rxn_mapping[reaction.id]) == 1:
                gene_asso = rxn_mapping[reaction.id]
            else:
                gene_asso = ['({})'.format(gene) for gene in rxn_mapping[reaction.id]]
            outfile.write('{}\t{}\n'.format(reaction.id, ' or '.join(gene_asso)))

    # Create a model_def file with all of the new reactions
    with open('./{}/model_def.tsv'.format(out), mode='w') as f:
        for reaction in out_nm.reactions:
            f.write('{}\n'.format(reaction.id))

    # Create a reactions.yaml file
    with open('./{}/reactions.yaml'.format(out), 'w') as out:
        out_mw.write_reactions(out, out_nm.reactions)

    # Create a compounds.yaml file
    with open('./{}/compounds.yaml'.format(out), 'w') as out:
        out_mw.write_reactions(out, out_nm.compounds)

    # Create a model.yaml file
    with open('./{}/model.yaml'.format(out), mode='w') as myam:
        myam.write('default_flux_limit: 100\n')
        myam.write('default_compartment: {}\n'.format("c"))
        myam.write('extracellular: null\n')
        myam.write('biomass: null\n')
        myam.write('compounds:\n')
        myam.write('- include: ./compounds.yaml\n')
        myam.write('reactions:\n')
        myam.write('- include: ./reactions.yaml\n')
        myam.write('- include: ./gene-association.tsv\n')
        myam.write('  format: tsv\n')
        myam.write('model:\n')
        myam.write('- include: model_def.tsv\n')

