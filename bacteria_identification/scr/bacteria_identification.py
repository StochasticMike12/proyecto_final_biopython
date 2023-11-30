'''
NAME

    bacteria_identification.py 

VERSION

    1.0.1

AUTHORES

    Camila Villazón Soto Innes
    Miguel Ángel Flores Varela

CONTACT

    mfvarela@lcg.unam.mx
    camilav@lcg.unam.mx

DESCRIPTION

    Este programa permite identificar el gnénero de bacterias de acuerdo 
    a los datos y parámetros establecidos por el usuario. El programa puede
    recibir como input un archivo fasta con la secuencia del gen 16S de la(s)
    bacteria(s) a identificar o bien un archivo xml con los resultados de 
    haber corrido BLAST con la(s) secuencia(s), esto último con la intención
    de agilizar el proceso de identificación.

CATEGORY

    Genómica

USAGE

    % py baacteria_identification.py <path_archivo_de_entrada> <fomarto>
  
'''


# ===========================================================================
# =                            imports
# ===========================================================================

# Importar librerías necesarias.

from Bio.Blast import NCBIXML
import identifier as ident
from Bio import SeqIO
import argparse


# ===========================================================================
# =                            Command Line Options
# ===========================================================================

# Definir los argumentos.

parser= argparse.ArgumentParser(description='Este programa permite identificar el gnénero de bacterias de acuerdo a los datos y parámetros establecidos por el usuario. El programa recibe como primer argumento el path del archivo a analizar y se especifica su formato (fasta o xml) como segundo argumento.')

parser.add_argument('Path',
                    metavar='path',
                    type=str,
                    help='Path del archivo')

parser.add_argument('Format',
                    metavar='format',
                    type=str,
                    help='Formato del archivo')

# Ejecutar método parse_args()
args = parser.parse_args()



# ===========================================================================
# =                            main
# ===========================================================================

# Abrir archivo de entrada. 
try:
    # En caso de que el archivo de entrada sea fasta, se ejecuta funcion
    # makeblast para correr BLAST
    if args.Format=='fasta':
        for record in SeqIO.parse(args.Path, "fasta"):
            ident.makeblast(record)
            file=record.id+'.xml'
            # Correr analisis de datos despues de correr BLAST
            n=1
            for blast_record in NCBIXML.parse(open(file)):
                ident.identify(blast_record, n)
                n+=1
    # Correr análisis de datos de archivos xml
    elif args.Format=='xml':
        n=1
        for record in NCBIXML.parse(open(args.Path)):
            ident.identify(record, n)
            n+=1
# Mensaje de error en caso de no encontrar archivo
except FileNotFoundError:
    print("No se encontró el archivo ingresado.")
