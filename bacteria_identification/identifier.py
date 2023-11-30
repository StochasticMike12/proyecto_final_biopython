from Bio.Blast import NCBIWWW
from Bio import Entrez, SeqIO 
import pandas as pd
import re



def makeblast(record):
    NCBIWWW.email='mfvarela@lcg.unam.mx'
    result_handle=NCBIWWW.qblast("blastn", 
                                 "nt",
                                 record.format("fasta"), 
                                 url_base='https://blast.ncbi.nlm.nih.gov/Blast.cgi',
                                 alignments=10,
                                 expect=0.05,
                                 format_object='Alignment',
                                 )
    file=record.id+'.xml'
    with open(file,"w") as save_to:
        save_to.write(result_handle.read())
        result_handle.close()

def identify(record):
    Entrez.email = "mfvarela@lcg.unam.mx"
    pidents=[]
    accessions=[]

    for alignment in record.alignments:
        pident=(alignment.hsps[0].identities/alignment.hsps[0].align_length)*100
        pidents.append(pident)
        accessions.append(alignment.accession)

    hits = pd.DataFrame.from_dict({"pident": pidents, "accession": accessions}).sort_values(by="pident")
    hits=hits[::-1]

    pidents=list(hits['pident'])

    if pidents[0]<99:
        print("\n\n> No hay ningún hit con porcentage de identidad mayor a 99%, por lo que es probable que se trate de una nueva especie.")
    
    print("\n\n> Los 5 mejores hits de acuerdo al porcentaje de identidad corresponden a los siguientes organismos:\n\n")
    
    accessions=list(hits['accession'])[0:5]
    pidents=pidents[0:5]
    organisms=[]

    for accesion in accessions:
        handle = Entrez.efetch(db="nucleotide", id=accesion, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        organisms.append(record.annotations['organism'])

    best_hits = pd.DataFrame.from_dict({"Organismo": organisms, "Porcetaje de identidad": pidents})
    print(best_hits)
    print("\n\n> Géneros a los que puede corresponder tu secuencia:\n")
    genera=[]
    for organism in best_hits['Organismo']:
        genera.append(re.findall(r'^[A-Z]\w+',organism)[0])
    genera=list(set(genera))
    for genus in genera:
        print("\n\t-",genus)
    print('\n')


