"""Este módulo permite la identificación del género de bacterias mediante 
funciones que blastean y analizan los datos ingreados."""

# ===========================================================================
# =                            imports
# ===========================================================================

# Importar librerías necesarias.

from Bio.Blast import NCBIWWW
from Bio import Entrez, SeqIO 
import pandas as pd
import re
import matplotlib.pyplot as plt



# ===========================================================================
# =                            functions
# ===========================================================================


def makeblast(record):
    """Esta función realiza un blast de la secuencia que recibe"""
    
    # Acceder a NCBI con un correo
    NCBIWWW.email='mfvarela@lcg.unam.mx'

    # Correr BLAST con los parámetros establecidos
    result_handle=NCBIWWW.qblast("blastn", 
                                 "nt",
                                 record.format("fasta"), 
                                 url_base='https://blast.ncbi.nlm.nih.gov/Blast.cgi',
                                 alignments=10,
                                 expect=0.05,
                                 format_object='Alignment',
                                 )
    
    # Establecer el nombre del xml
    file=record.id+'.xml'

    # Guardar los resultados del BLAST en el archivo xml
    with open(file,"w") as save_to:
        save_to.write(result_handle.read())
        result_handle.close()


def identify(record, n):
    """Esta función analiza los datos obtenidos del output de BLAST para 
    identificar el género de las bacterias de las cuales se ingresó la 
    secuencia de su gen 16S."""
    
    # Acceder a entrez con un correo
    Entrez.email = "mfvarela@lcg.unam.mx"
    
    # Imprimir el número de secuencia y su id
    print("\n\n> Secuencia", n,":",record.query)

    # Obtener el porcentaje de identidad y el id de cada hit
    pidents=[]
    accessions=[]
    for alignment in record.alignments:
        pident=(alignment.hsps[0].identities/alignment.hsps[0].align_length)*100
        pidents.append(pident)
        accessions.append(alignment.accession)

    # Ordenar los hits de acuerdo al porcentaje de identidad
    hits = pd.DataFrame.from_dict({"pident": pidents, "accession": accessions}).sort_values(by="pident")
    hits=hits[::-1]

    # Analizar si existen posibles nuevas especies de cuerdo al porcentaje de identidad
    pidents=list(hits['pident'])
    if pidents[0]<99:
        print("\n\n> No hay ningún hit con porcentage de identidad mayor a 99%, por lo que es probable que se trate de una nueva especie.")
    
    # Obtener los códigos de acceso de los mejores 5 hits de acuerdo al porcentaje de identidad
    print("\n\n> Los 5 mejores hits de acuerdo al porcentaje de identidad corresponden a los siguientes organismos:\n\n")
    accessions=list(hits['accession'])[0:5]
    pidents=pidents[0:5]
    
    # Consultar a qué organismo pertenecen los mejores 5 hits
    organisms=[]
    for accesion in accessions:
        handle = Entrez.efetch(db="nucleotide", id=accesion, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        organisms.append(record.annotations['organism'])

    # Imprimir y graficar el organismos al correspondiente a cada hit junto con su respectivo porcentaje de identidad
    bacteria=query(organisms, pidents)
    bacteria.get_best_hits()
    print(bacteria.best_hits)
    bacteria.plot()

    # Obtener los géneros de los mejores 5 hits
    print("\n\n> Géneros a los que puede corresponder tu secuencia:\n")
    genera=[]
    for organism in bacteria.best_hits['Organismo']:
        genera.append(re.findall(r'^[A-Z]\w+',organism)[0])
    genera=list(set(genera))
    for genus in genera:
        print("\n\t-",genus)
    print('\n')



# ===========================================================================
# =                            main
# ===========================================================================


class query():
    """Esta clase posee atributos y funciones que permiten visualizar los
    mejores hits de acuerdo al porcentaje de identidad."""
    
    def __init__(self, organims, pidents):
        """Atributos de la clase, es la información correspondiente a los 
        mejores 5 hits de acuerdo al porcentaje de identidad."""
        self.organims=organims
        self.pidents=pidents
        self.best_hits=None

    def get_best_hits(self):
        """Esta función crea un data frame con los el nombre de los 
        organismos y sus porcentajes de identidad correspondientes a los 
        mejores 5 hit"""
        self.best_hits=pd.DataFrame.from_dict({"Organismo": self.organims, "Porcentaje de identidad": self.pidents})

    def plot(self):
        """Esta función permite plotear los organismos de los mejores 5 
        hits de acuerdo a su porcentaje de identidad."""
        x=self.best_hits['Organismo']
        y=self.best_hits['Porcentaje de identidad']/100
        fig, ax = plt.subplots()
        ax.bar(x = x, height = y)
        plt.xlabel('Organismos')
        plt.ylabel('Porcentaje de identidad / 100')
        plt.show()
