import wget
import os
import re
import gzip
from random import randint
from collections import defaultdict
import subprocess
from datetime import datetime
import shutil

def download_references(accession, out):
    # set up the output directories
    os.mkdir(out)
    os.mkdir("".join([out,"/NCBI_references"]))
    
    # Use the accession number to build a list of URLs to downlaod from
    listall=re.split("_", accession)
    # Check the accession is valid for genbank or refseq
    if listall[0]!="GCA" and listall[0]!="GCF":
        raise Exception('This does not appear to be a valid accession number!')
    url="/".join(["ftp://ftp.ncbi.nlm.nih.gov/genomes/all",listall[0],listall[1][0:3],
        listall[1][3:6],listall[1][6:9],accession,accession])
    proteome="".join([url,"_protein.faa.gz"])
    genome="".join([url,"_genomic.fna.gz"])
    print(genome)
    gff="".join([url,"_genomic.gff.gz"])
    try:
        genome_out = wget.download(genome, out="".join([out,"/NCBI_references"]))
    except:
        raise Exception('Genome Download Failed! Are we sure this is a valid accession?')
    try:
        proteome_out = wget.download(proteome, out="".join([out,"/NCBI_references"]))
    except:
        raise Exception('Proteome Download Failed! Are we sure this is a valid accession?')
    try:
        gff_out = wget.download(gff, out="".join([out,"/NCBI_references"]))
    except:
        raise Exception('GFF Download Failed! Are we sure this is a valid accession?')
    genome="".join([out,"/NCBI_references/",accession,"_genomic.fna.gz"])
    proteome="".join([out,"/NCBI_references/",accession,"_protein.faa.gz"])
    gff="".join([out,"/NCBI_references/",accession,"_genomic.gff.gz"])
    return(genome, proteome, gff)

def run_diamond(assembly, proteome, out):
    # Create a directory to store all of the diamond files
    diam_out="".join([out,"/diamond_out/"])
    os.mkdir(diam_out)
    
    # Unzip the proteome
    logfile=open("".join([diam_out,"logfile.txt"]),"a")
    #run_gunzip = subprocess.run(['gunzip -f ', proteome], shell=True, stdout=logfile, stderr=subprocess.STDOUT)
    print(proteome)
    with gzip.open(proteome, 'rb') as infile:
        with open(proteome[:-3], 'wb') as outfile:
            shutil.copyfileobj(infile, outfile)
    proteome=proteome[:-3]
    print(proteome)

    # make the diamond database from the reference
    db=''.join([diam_out, "diamond_db"])
    print('diamond makedb --in ',proteome,' --db ',db)
    run_makedb = subprocess.run("".join(['diamond makedb --in ',proteome,' --db ',db]), shell=True, stdout=logfile, stderr=subprocess.STDOUT)
    db=''.join([diam_out, "diamond_db.dmnd"])
    output=''.join([diam_out, "diamond.tsv"])
    run_diamond = subprocess.run("".join(['diamond blastp --query ',assembly,' --db ',db,' --out ',output,' --outfmt 6 qseqid sseqid pident evalue length qlen slen qstart qend sstart send sseq']), shell=True, stdout=logfile, stderr=subprocess.STDOUT)
    logfile.close()
    return(output)

def create_prot(accession, out, assembly, pident, pcover):
    # Time it! Why not?
    startTime=datetime.now()

    # Download the reference genome, proteome, and gff files.
    genome, proteome, gff = download_references(accession, out)
    
    # Loop through the genome file and write the first lines of 
    # the karyotype file, which denote the contig lengths and colors.
    # Also outputs randomly generated colors for each contig based on
    # rbg(0,0,0) format.
    key='DONT'
    sequence=''
    label=0
    color_dict={}
    with gzip.open(genome, "rb") as infile:
        with open("".join([out,"/",accession, ".txt"]), "w") as outfile:
            for line in infile:
                line=str(line, 'utf-8')
                line=line.rstrip()
                if line[0]==">":
                    # This is included simply to skip the first iteration
                    if key != "DONT":
                            outfile.write(" ".join(["chr -",key, str(label), "0", str(len(sequence)), "".join([key, "_col"])]))
                            outfile.write("\n")
                            color_dict["".join([key, "_col"])] = ",".join([str(r),str(b),str(g)])
                    label+=1
                    r=randint(0,255)
                    b=randint(0,255)
                    g=randint(0,255)
                    listall=re.split("\s", line)
                    key=listall[0][1:]
                    sequence=''
                else:
                    sequence+=line.strip()
            outfile.write(" ".join(["chr -",key, str(label), "0", str(len(sequence)), "".join([key, "_col"])]))
            outfile.write("\n")
            color_dict["".join([key, "_col"])] = ",".join([str(r),str(b),str(g)])
    # Export the color coding
    with open("".join([out,"/",accession, "_colors.conf"]), "w") as outfile:
        for i in color_dict:
            outfile.write("\t".join([i, "=", color_dict[i]]))
            outfile.write("\n")
    
    # Loop through the gff file and find the contig and star/end
    # position for every gene. Stores as a dictionary for output
    # after cross referenceing with the proteome file. Make the
    # genes one color
    r=randint(0,255)
    b=randint(0,255)
    g=randint(0,255)
    color_proteome=",".join([str(r),str(b),str(g)])
    band_proteome={}
    position_dict={}
    protein_total=0
    with gzip.open(gff, "rb") as infile:
        for line in infile:
            line=str(line, 'utf-8')
            if line[0]!="#":
                listall_space=re.split("\t", line)
                listall_semi=re.split(";", listall_space[8])
                if listall_semi[0][0:7]=="ID=cds-":
                    band_proteome[listall_semi[0][7:]]=" ".join([listall_space[0],listall_space[3],listall_space[4]])
                    position_dict[listall_semi[0][7:]]=(listall_space[0],listall_space[3],listall_space[4])
    with gzip.open(proteome, "rb") as infile:
        with open("".join([out,"/",accession, "_bacteriaFeature.txt"]), "w") as outfile:
            for line in infile:
                line=str(line, 'utf-8')
                if line[0]==">":
                    listall=re.split("\s", line)
                    protein_total+=1
                    try:
                        outfile.write(band_proteome[listall[0][1:]])
                        outfile.write("\n")
                    except:
                        sys.stdout.write("".join(["""your code sucks, you missed a protein (or more likely you corrupted one of the
                        input files somehow): """, listall[0]]))
                        exit(0)
    with open("".join([out,"/",accession, "_colors.conf"]), "w") as outfile:
        outfile.write("\t".join(["color_proteome", "=", color_proteome]))
        outfile.write("\n")

    # Run diamond on the assembly
    diam_out = run_diamond(assembly, proteome, out)

    # Use the diamond output to generate the foram proteome.
    # Parse out the blast table and generate a highlight for all of
    # the features that meet the specified threshold. Then generate
    # a feature file.
    protein_mapped=set()
    proteome_id=defaultdict(lambda:0)
    with open(diam_out, "r") as infile:
        with open("".join([out,"/",accession, "_foramFeature.txt"]), "w") as outfile:
            for line in infile:
                listall=re.split("\t", line)
                pcover_calc=(float(listall[4])) / float(listall[6]) * 100
                if float(listall[2]) >= int(pident) and pcover_calc >= int(pcover):
                    try:
                        contig, start, end = position_dict[listall[1]]
                        outfile.write(" ".join([contig, start, end]))
                        outfile.write("\n")
                        protein_mapped.add(listall[1])
                        # extract the best hits
                        if float(listall[2]) > float(proteome_id[listall[1][1]]):
                            proteome_id[listall[1]]=[listall[0],float(listall[2])]
                    except:
                        pass

    # Use the protein_mapped variable to generate a list of the 
    # unique bacterial proteins that were mapped to
    with open("".join([out,"/",accession, "_bacterial_names.tsv"]), "w") as outfile:
        for i in protein_mapped:
            outfile.write(i)
            outfile.write("\n")

    # Use the proteome_id variable to generate a fasta file from 
    # the assembly.
    # First generate the reference database from the assembly.
    key=''
    sequence=''
    seqs={}
    with open(assembly, 'r') as seq_file:
        for line in seq_file:
            line = line.strip()
            if line.startswith('>'):
                sequence=sequence.rstrip()
                seqs[key]=sequence
                listall=re.split("\s", line)
                key=listall[0][1:]
                sequence=''
            else:
                sequence+=line.strip()
    # Generate the fasta file based on the contig names present
    # in the proteome_id variable
    with open("".join([out,"/",accession, "_foram_proteome.fa"]), "w") as outfile:
        for i in proteome_id:
            if proteome_id[i]!=0:
                sequence=seqs[proteome_id[i][0]]
                outfile.write(">" + proteome_id[i][0] + "\n" + sequence + "\n")

    # Generate the circos configuration file.
    with open("".join([out,"/",accession, ".conf"]), "w") as outfile:
        outfile.write("\t".join(["karyotype", "=", "".join([out, ".txt\n"])]))
        outfile.write("\n".join(["<image>"," radius = 1500p"," background   = white\n"]))
        outfile.write("\t".join([" file*    = ","".join([out, "_foram.png\n"])]))
        outfile.write("</image>\n")
        outfile.write("\n") 
        outfile.write("\n".join(["<highlights>"," z = 5"," <highlight>\n"]))
        outfile.write("\t".join([" file","=","".join([out, "_foramFeature.txt\n"])]))
        outfile.write("\n".join([" fill_color       = lred"," r0    = 0.5r"," r1    = 0.8r", " </highlight>\n"]))
        outfile.write("\n")
        outfile.write(" <highlight>\n")
        outfile.write("\t".join([" file","=","".join([out, "_bacteriaFeature.txt\n"])]))
        outfile.write("\n".join([" fill_color       = lblue"," r0    = 0.8r"," r1    = 1.0r", " </highlight>","</highlights>\n"]))
        outfile.write("\n")
        outfile.write("\n".join(["<ideogram>"," show_bands  = yes"," fill_bands     = yes"," band_transparency      = 4",
            " band_stroke_thickness = 2"," show_label       = yes"," label_radius     = dims(ideogram,radius) + 0.075r", 
            " label_size       = 24"," label_parallel   = yes\n"," <spacing>"," default     = 0.005r"," </spacing>\n\n",
            " radius        = 0.90r"," thickness        = 200p"," fill             = yes"," stroke_color     = dgrey",
            " stroke_thickness = 2p","</ideogram>\n\n"]))
        outfile.write("\n".join(["<image>"," <<include etc/image.conf>>","</image>","<<include etc/colors_fonts_patterns.conf>>\n"]))
        outfile.write("".join(["<<include ","".join([out, "_colors.conf"]),">>\n"]))
        outfile.write("<<include etc/housekeeping.conf>>")

    # Generate run statistics
    with open("".join([out,"/","Run_Statistics.txt"]), "w") as outfile:
        outfile.write("".join(["# Runtime: ", str(datetime.now() - startTime), "\n\n"]))
        outfile.write("# The data in this directory was parsed at: \n")
        outfile.write("".join(["# ",str(pident), "% ID\n"]))
        outfile.write("".join(["# ",str(pcover), "% Cover\n\n"]))
        outfile.write(" ".join(["# The total number of proteins in the proteome is:", str(protein_total), "\n"]))
        outfile.write(" ".join(["# The number of proteins mapped to from the assembly is:", str(len(protein_mapped)), "\n\n"]))
        outfile.write("# For each protein in the foram proteome, this is the % ID match:\n")
        for i in proteome_id:
            if proteome_id[i]!=0:
                outfile.write("\t".join([proteome_id[i][0],i,str(proteome_id[i][1])]))
                outfile.write("\n")
