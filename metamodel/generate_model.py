import sys
import re
import os
from psamm.datasource.native import NativeModel, ModelReader, ModelWriter
from psamm.expression import boolean
from collections import defaultdict
from psamm.datasource.context import FileMark
from psamm.datasource.reaction import Reaction, Compound
from psamm.expression.affine import Expression
from psamm.formula import Formula
import yaml
from collections import OrderedDict
from psamm.datasource.native import parse_reaction_equation_string
from Bio.KEGG import REST
from Bio.KEGG import Enzyme

def dict_representer(dumper, data):
    return dumper.represent_dict(data.items())

def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

'''
This function will create a draft model that is essentially
a subset of an existing model based on the reactions provided
in the reactions:genes file input. Typically, this existing
model we use is one based on all of the reactions in kegg,
but realistically you could provide this with any yaml
formatted model.yaml
'''
def create_model_kegg(out, rxn, kegg):
    os.mkdir(out)

    # Initialize some stuff that store the reaction data
    # Stores the Reaction:Gene associations
    rxn_mapping=defaultdict(list)
    with open(rxn, "r") as infile:
        for line in infile:
            line=line.rstrip()
            listall=re.split("\t", line)
            if len(listall) > 1:
                genes=re.split(",", listall[1])
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

    # Check through compounds in the kegg model. If a needed compound
    # add to the new model.
    for compound in nm.compounds:
        if compound.id in compound_list:
            out_nm.compounds.add_entry(compound)

    # Create a model writer object that will write the new native
    # model out
    out_mw = ModelWriter()

    with open('{}/gene-association.tsv'.format(out), mode='w') as outfile:
        outfile.write('id\tgenes\n')
        for reaction in out_nm.reactions:
            if len(rxn_mapping[reaction.id]) == 1:
                gene_asso = rxn_mapping[reaction.id]
            else:
                gene_asso = ['({})'.format(gene) for gene in rxn_mapping[reaction.id]]
            outfile.write('{}\t{}\n'.format(reaction.id, ' or '.join(gene_asso)))

    # Create a model_def file with all of the new reactions
    with open('{}/model_def.tsv'.format(out), mode='w') as f:
        for reaction in out_nm.reactions:
            f.write('{}\n'.format(reaction.id))

    # Create a reactions.yaml file
    with open('{}/reactions.yaml'.format(out), 'w') as f:
        out_mw.write_reactions(f, out_nm.reactions)

    # Create a compounds.yaml file
    with open('{}/compounds.yaml'.format(out), 'w') as f:
        out_mw.write_reactions(f, out_nm.compounds)

    # Create a model.yaml file
    with open('{}/model.yaml'.format(out), mode='w') as myam:
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

'''
This function creates a draft model based on the reactions:genes
file specified in the rxn variable. This function generates
reaction and compoiund information by utilizing the kegg
REST api to download the reaction information and uses psamm
functions to parse out the kegg data. This is the default
function, as not everyone will have a reference model to base
their current model on.
'''
def create_model_api(out, rxn):
    os.mkdir(out)

    # Initialize some stuff that store the reaction data
    # Stores the Reaction:Gene associations
    rxn_mapping=defaultdict(list)
    with open(rxn, "r") as infile:
        for line in infile:
            line=line.rstrip()
            listall=re.split("\t", line)
            if len(listall) > 1 and line[0]!="G":
                genes=re.split(",", listall[1])
                rxn_mapping[listall[0]]=genes

    with open(os.path.join(out, "log.tsv"), "a+") as f:
        f.write("List of invalid Kegg IDs: \n")
    # Generate the yaml file for the reactions
    reaction_entry_list = []
    for entry in _parse_kegg_entries(out, rxn_mapping, ReactionEntry):
        reaction_entry_list.append(entry)

    # Sets up the yaml object for the reactions and writes
    # out the parsed reaction information to a reactions.yaml
    # file
    yaml.add_representer(OrderedDict, dict_representer)
    yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
                         dict_constructor)
    yaml_args = {'default_flow_style': False,
                 'encoding': 'utf-8',
                 'allow_unicode': True}

    # Generate the yaml file for the compounds
    compound_entry_list=[]
    generic_entry_list=[]
    compound_set=set()
    global generic_compounds_list
    generic_compounds_list=[]
    for reaction in reaction_entry_list:
        eq = parse_reaction_equation_string(reaction.equation, 'c')
        for i in eq.compounds:
            compound_set.add(str(i[0].name))
    for entry in _parse_kegg_entries(out, compound_set, CompoundEntry):
        compound_entry_list.append(entry)
    with open(os.path.join(out, 'compounds.yaml'), 'w+') as f:
        yaml.dump(list(model_compounds(compound_entry_list)), f, **yaml_args)

    for entry in _parse_kegg_entries(out, generic_compounds_list, CompoundEntry):
        generic_entry_list.append(entry)
    with open(os.path.join(out, 'compounds_generic.yaml'), 'w+') as f:
        yaml.dump(list(model_generic_compounds(generic_entry_list)), f, **yaml_args)
    with open(os.path.join(out, 'log.tsv'), "a+") as f:
        f.write("\nThere are {} generic compounds in the model\n".format(str(len(generic_compounds_list))))
        f.write("Generic compounds:\n")
        for i in generic_entry_list:
            f.write("{}".format(i.id))
            if i.name:
                f.write("|{}".format(i.name))
            else:
                f.write("|")
            if i.formula:
                f.write("|{}\n".format(i.formula))
            else:
                f.write("|\n")
    reaction_list_out = []
    reaction_list_generic = []
    with open(os.path.join(out, 'log.tsv'), 'a+') as f:
        f.write("\nThe reactions containing these generic compounds are: \n")
        for reaction in reaction_entry_list:
            if any(i in str(reaction.equation) for i in generic_compounds_list):
                f.write("{}".format(reaction.id))
                if reaction.name:
                    f.write("|{}".format(reaction.name))
                else:
                    f.write("|")
                if reaction.equation:
                    f.write("|{}\n".format(reaction.equation))
                else:
                    f.write("|\n")
                reaction_list_generic.append(reaction)
            else:
                reaction_list_out.append(reaction)
    with open(os.path.join(out, 'reactions.yaml'), 'w+') as f:
        yaml.dump(list(model_reactions(reaction_list_out)), f, **yaml_args)
    with open(os.path.join(out, 'reactions_generic.yaml'), 'w+') as f:
       	yaml.dump(list(model_reactions(reaction_list_generic)), f, **yaml_args)


    # Create a model.yaml file
    with open('{}/model.yaml'.format(out), mode='w') as myam:
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

    with open('{}/gene-association.tsv'.format(out), mode='w') as outfile:
        outfile.write('id\tgenes\n')
        for reaction in reaction_list_out:
            if len(rxn_mapping[reaction.id]) == 1:
                gene_asso = rxn_mapping[reaction.id]
            else:
                gene_asso = ['({})'.format(gene) for gene in rxn_mapping[reaction.id]]
            outfile.write('{}\t{}\n'.format(reaction.id, ' or '.join(gene_asso)))

    with open('{}/gene-association_generic.tsv'.format(out), mode='w') as outfile:
        outfile.write('id\tgenes\n')
        for reaction in reaction_list_generic:
            if len(rxn_mapping[reaction.id]) == 1:
                gene_asso = rxn_mapping[reaction.id]
            else:
                gene_asso = ['({})'.format(gene) for gene in rxn_mapping[reaction.id]]
            outfile.write('{}\t{}\n'.format(reaction.id, ' or '.join(gene_asso)))

    # Create a model_def file with all of the new reactions
    with open('{}/model_def.tsv'.format(out), mode='w') as f:
        for reaction in reaction_list_out:
            f.write('{}\n'.format(reaction.id))

    # Write out some final statistics for the model
    with open(os.path.join(out, 'log.tsv'), 'a+') as f:
        f.write("\nThere are {} reactions in the model".format(str(len(reaction_list_out))))
       	f.write("\nThere are {} compounds in the model\n".format(str(len(compound_entry_list))))

'''
Function to sort the downloaded kegg object into a format
that is compatible with the psamm api for storage in
a reactions.yaml file.
'''
def model_reactions(reaction_entry_list):
        for reaction in reaction_entry_list:
            d = OrderedDict()
            d['id'] = encode_utf8(reaction.id)

            enzymes_list = []
            for i in reaction.enzymes:
                enzymes_list.append(i)

            pathways_list = []
            if reaction.pathways is not None:
                for i in reaction.pathways:
                    pathways_list.append(i[1])
            if len(pathways_list) == 0:
                pathways_list = None

            orth_list = []
            for i in reaction.orthology:
                    orth_list.append(i)
            if len(orth_list) == 0:
                orth_list = None

            rpairs_list = []
            for i in reaction.rpairs:
                rpairs_list.append(i)

            if hasattr(reaction, 'name') and reaction.name is not None:
                d['name'] = encode_utf8(reaction.name)
            if hasattr(reaction, 'names') and reaction.names is not None:
                names_l = []
                for i in reaction.names:
                    names_l.append(i)
                d['names'] = encode_utf8(names_l)
            if hasattr(reaction, 'equation') and reaction.equation is not None:
                d['equation'] = encode_utf8(str(reaction.equation))
            if hasattr(reaction, 'enzymes') and reaction.enzymes is not None:
                d['enzymes'] = encode_utf8(enzymes_list)
            if hasattr(reaction, 'pathways') and reaction.pathways is not None:
                d['pathways'] = encode_utf8_list(pathways_list)
            if hasattr(reaction, 'comment') and reaction.comment is not None:
                d['comment'] = encode_utf8(str(reaction.comment))
            if hasattr(reaction, 'orthology') and reaction.orthology is not None:
                d['orthology'] = encode_utf8(orth_list)
            if hasattr(reaction, 'rpairs') and reaction.rpairs is not None:
                d['rpairs'] = encode_utf8(rpairs_list)
            yield d

'''
Function to sort the downloaded	kegg object into a format
that is	compatible with the psamm api for storage in
a compounds.yaml file.
'''
def model_compounds(compound_entry_list):
        non_gen_compounds = []
        global generic_compounds_list
        for compound in compound_entry_list:
            try:
                form = Formula.parse(str(compound.formula))
                if form.is_variable():
                    generic_compounds_list.append(compound.id)
                    continue
                elif compound.formula is None:
                    generic_compounds_list.append(compound.id)
                    continue
                elif 'R' in str(compound.formula):
                    generic_compounds_list.append(compound.id)
                    continue
                else:
                    d = OrderedDict()
                    d['id'] = encode_utf8(compound.id)
                    non_gen_compounds.append(compound.id)
                    if hasattr(compound, 'name') and compound.name is not None:
                        d['name'] = encode_utf8(compound.name)
                    if hasattr(compound, 'names') and compound.names is not None:
                        names_l = []
                        for i in compound.names:
                            names_l.append(i)
                        d['names'] = encode_utf8(names_l)
                    if hasattr(compound, 'formula') and compound.formula is not None:
                        d['formula'] = encode_utf8(str(compound.formula))
                    if hasattr(compound, 'mol_weight') and compound.mol_weight is not None:
                        d['mol_weight'] = encode_utf8(compound.mol_weight)
                    if hasattr(compound, 'comment') and compound.comment is not None:
                        d['comment'] = encode_utf8(str(compound.comment))
                    if hasattr(compound, 'dblinks') and compound.dblinks is not None:
                        for key, value in compound.dblinks:
                            d['{}'.format(key)] = encode_utf8(value)
                    yield d
            except:
                generic_compounds_list.append(compound.id)

def model_generic_compounds(compound_entry_list):
        #generic_compounds_list = []
        #global generic_compounds_list
        non_gen_compounds = []
        for compound in compound_entry_list:
            try:
                form = Formula.parse(str(compound.formula))
                d = OrderedDict()
                d['id'] = encode_utf8(compound.id)
                non_gen_compounds.append(compound.id)
                if hasattr(compound, 'name') and compound.name is not None:
                    d['name'] = encode_utf8(compound.name)
                if hasattr(compound, 'names') and compound.names is not None:
                    names_l = []
                    for i in compound.names:
                        names_l.append(i)
                    d['names'] = encode_utf8(names_l)
                if hasattr(compound, 'formula') and compound.formula is not None:
                    d['formula'] = encode_utf8(str(compound.formula))
                if hasattr(compound, 'mol_weight') and compound.mol_weight is not None:
                    d['mol_weight'] = encode_utf8(compound.mol_weight)
                if hasattr(compound, 'comment') and compound.comment is not None:
                    d['comment'] = encode_utf8(str(compound.comment))
                if hasattr(compound, 'dblinks') and compound.dblinks is not None:
                    for key, value in compound.dblinks:
                        d['{}'.format(key)] = encode_utf8(value)
                yield d
            except:
                print(compound)
                continue


'''
Two functions listed here for compatibility with
python 2 and python 3
'''
def encode_utf8(s):
    is_python2 = sys.version_info.major == 2
    if is_python2:
        if isinstance(s, unicode):
            return s.encode('utf-8')
    else:
        return s
def encode_utf8_list(s):
    is_python2 = sys.version_info.major == 2
    if is_python2:
        if isinstance(s, unicode):
            return s[1].encode('utf-8')
    else:
        return s

'''
Downloads the kegg entry associated with a reaction or
a compound and stores each line in an object that can
be parsed as a reaction or a compound, depending on the
input
'''
def _parse_kegg_entries(out, rxn_mapping, entry_class, context=None):
    # go through the rxn mapping dict and pull all of the
    # reaction information out of the Kegg API
    with open(os.path.join(out,'log.tsv'), "a+") as f:
        for reactions in rxn_mapping:
            if reactions[0]=="R" or reactions[0]=="C":
                try:
                    request = REST.kegg_get(reactions)
                except:
                    f.write("".join(["  - ",reactions,"\n"]))
                    continue
                entry_line = None
                section_id = None
                reaction = {}
                for lineno, line in enumerate(request):
                    line=line.rstrip()
                    if line=="///":
                        continue
                    if entry_line is None:
                        entry_line = lineno
                    # Look for the beginning of the section
                    m = re.match(r'([A-Z_]+)\s+(.*)', line.rstrip())
                    if m is not None:
                        section_id = m.group(1).lower()
                        reaction[section_id] = [m.group(2)]
                    elif section_id is not None:
                        reaction[section_id].append(line.strip())
                    else:
                        raise ParseError(
                            'Missing section identifier at line {}'.format(lineno))
                mark = FileMark(context, entry_line, 0)
                yield entry_class(reaction, filemark=mark)

class ReactionEntry(object):
    """Representation of entry in KEGG reaction file"""

    def __init__(self, values, filemark=None):
        self.values = dict(values)
        if 'entry' not in values:
            raise ParseError('Missing reaction identifier')
        self._id, _ = values['entry'][0].split(None, 1)
        self._filemark = filemark

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        try:
            return next(self.names)
        except StopIteration:
            return None

    @property
    def names(self):
        if 'name' in self.values:
            for line in self.values['name']:
                for name in line.rstrip(';').split(';'):
                    yield name.strip()

    @property
    def definition(self):
        if 'definition' not in self.values:
            return None
        return self.values['definition'][0]

    @property
    def equation(self):
        if 'equation' not in self.values:
            return None
        return self.values['equation'][0]

    @property
    def enzymes(self):
        if 'enzyme' in self.values:
            for line in self.values['enzyme']:
                for enzyme in line.split():
                    yield enzyme

    @property
    def pathways(self):
        if 'pathway' in self.values:
            for line in self.values['pathway']:
                pathway, name = line.split(None, 1)
                yield pathway, name

    @property
    def comment(self):
        if 'comment' not in self.values:
            return None
        return '\n'.join(self.values['comment'])

    @property
    def rpairs(self):
        #class parsing needs to be fixed at some point.
        if 'rclass' in self.values:
            r_dict = {}
            for line in self.values['rclass']:
                en = line.split(None, 2)
                pair, compound_l = en[0], en[1:]
                for compounds in compound_l:
                    compounds = list(compounds.split('_', 1))
                    r_dict[pair] = compounds
                    yield r_dict

    @property
    def orthology(self):
        if 'orthology' in self.values:
            KO_name_dict = {}
            for line in self.values['orthology']:
                split = line.split(None, 1)
                yield split[0]

    def __getitem__(self, name):
        if name not in self.values:
            raise AttributeError('Attribute does not exist: {}'.format(name))
        return self._values[name]

    @property
    def filemark(self):
        return self._filemark

    def __repr__(self):
        return '<ReactionEntry "{}">'.format(self.id)

class ParseError(Exception):
    """Exception used to signal errors while parsing"""

class CompoundEntry(object):
    """Representation of entry in KEGG compound file"""

    def __init__(self, values, filemark=None):
        self.values = dict(values)
        if 'entry' not in values:
            raise ParseError('Missing compound identifier')
        self._id, _ = values['entry'][0].split(None, 1)
        self._filemark = filemark

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        try:
            return next(self.names)
        except StopIteration:
            return None

    @property
    def names(self):
        if 'name' in self.values:
            for line in self.values['name']:
                for name in line.rstrip(';').split(';'):
                    yield name.strip()

    @property
    def reactions(self):
        if 'reaction' in self.values:
            for line in self.values['reaction']:
                for rxnid in line.split():
                    yield rxnid

    @property
    def enzymes(self):
        if 'enzyme' in self.values:
            for line in self.values['enzyme']:
                for enzyme in line.split():
                    yield enzyme

    @property
    def formula(self):
        if 'formula' not in self.values:
            return None
        return self.values['formula'][0]

    @property
    def exact_mass(self):
        if 'exact_mass' not in self.values:
            return None
        return float(self.values['exact_mass'][0])

    @property
    def mol_weight(self):
        if 'mol_weight' not in self.values:
            return None
        return float(self.values['mol_weight'][0])

    @property
    def pathways(self):
        if 'pathway' in self.values:
            for line in self.values['pathway']:
                pathway, name = line.split(None, 1)
                yield pathway, name

    @property
    def dblinks(self):
        if 'dblinks' in self.values:
            for line in self.values['dblinks']:
                database, entry = line.split(':', 1)
                yield database.strip(), entry.strip()

    @property
    def comment(self):
        if 'comment' not in self.values:
            return None
        return '\n'.join(self.values['comment'])

    def __getitem__(self, name):
        if name not in self.values:
            raise AttributeError('Attribute does not exist: {}'.format(name))
        return self._values[name]

    @property
    def filemark(self):
        return self._filemark

    def __repr__(self):
        return '<CompoundEntry "{}">'.format(self.id)
