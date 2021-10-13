import re
import sys
from collections import defaultdict

def parse_eggnog(out, orthology, type):
    # Dictionary of reactions to genes
    ortho_dict=defaultdict(list)
    # Populate the dictionary
    with open(orthology, "r") as infile:
        for line in infile:
            if line[0]!="#":
                line=line.rstrip()
                listall=re.split("\t", line)
                if type == "rxn":
                    if listall[11]!='':
                        rxns=re.split(",", listall[11])
                        for i in rxns:
                            ortho_dict[i].append(listall[0])
                elif type == "ec":
                    if listall[7]!='':
                        rxns=re.split(",", listall[7])
                        for i in rxns:
                            ortho_dict[i].append(listall[0])

    # generate the outfile
    with open(out, "w") as outfile:
        for i in ortho_dict:
            outfile.write("\t".join([i,",".join(ortho_dict[i])]))
            outfile.write("\n")
