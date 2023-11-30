# Bacteria Identification

Camila Villazón Soto Innes
Miguel Ángel Flores Varela

30/nov/2023

## Introducción
BLAST (Basic Local Aligment Search Tool) es una herramienta que permite buscar regiones de similitud entre secuencias biológicas, ya sen nucleotídicas o protéicas. Las secuencias de interés son comparadas con las bases de datos y se calcula la significancia estadística. Mediante un BLAST es posible inverir homologías o identificar miembros de familias de genes.

Entre estos datos estadísticos se destacan el porcentage de identidad (pident) que describe el porcentaje de matches idénticos en la secuencia; y el valor de expectancia (evalue) que nos habla de secuencias que se espera que coincidan por azar, es un valor que nos permite confiar en el match encontrado.

La maquinaria transcripcional de los organismos es altamente conservada. Una secuencia de especial importancia es la secuencia 16S de rRNA. con esta secuencia se pueden establecer relaciones filogenéticas o identificar a una bacteria. Por esta razón, la comparación de 16S de de gran importancia microbiológica. Cuando el pident sea mayor o igual al 99%, se considera la misma especie.

## Objetivo
Identificar genero de especies mediante 16S o posibles especies nuevas. 

## Método
1. Usuario ingresa el archivo de secuencias 16S en formato fasta o archivo xml con resultados de BLAST, indicando el formato del archivo.
2. Si el archivo está en formato fasta, el programa corre BLAST y el análisis de los resultados. En el caso contrario correrá únicamente el análisis de los resultados de BLAST.


## Resultados y pruebas

Caso 1. Archivo de entrada fasta
El programa se ejecuta `py bacteria_identification.py ..\data\pJET16S-1-4_pJET1_2F.fasta fasta`

El programa devuelve: un archivo xml con el outut de BLAST, plotea los mejores hits de acuerdo al porcentaje de identidad y los imprime junto con el género al que puede pertenecer la secuencia. En concreto el mensaje es el siguiente:
```
> Los 5 mejores hits de acuerdo al porcentaje de identidad corresponden a los siguientes organismos:


                       Organismo  Porcentaje de identidad
0            Bacillus velezensis                99.566349
1            Bacillus velezensis                99.566349
2  Bacillus sp. (in: firmicutes)                99.480519
3  Bacillus sp. (in: firmicutes)                99.480519
4              Bacillus subtilis                99.480519

> Géneros a los que puede corresponder tu secuencia:


        - Bacillus
```

Caso 2. Archivo de entrada xml
El programa se ejecuta `py bacteria_identification.py ..\data\pJET16S-1-4_pJET1_2.xml xml`

El programa devuelve: plotea los mejores hits de acuerdo al porcentaje de identidad y los imprime junto con el género al que puede pertenecer la secuencia. En concreto el mensaje es el siguiente:
```
> Secuencia 1 : H230705-030_O14_pJET16S-1-4_pJET1_2F.ab1

> Los 5 mejores hits de acuerdo al porcentaje de identidad corresponden a los siguientes organismos:


                       Organismo  Porcentaje de identidad
0            Bacillus velezensis                99.566349
1            Bacillus velezensis                99.566349
2  Bacillus sp. (in: firmicutes)                99.480519
3  Bacillus sp. (in: firmicutes)                99.480519
4              Bacillus subtilis                99.480519

> Géneros a los que puede corresponder tu secuencia:


        - Bacillus


> Secuencia 2 : H230705-030_M14_pJET16S-1-4_pJET1_2R.ab1

> Los 5 mejores hits de acuerdo al porcentaje de identidad corresponden a los siguientes organismos:


                    Organismo  Porcentaje de identidad
0         Bacillus velezensis                99.176277
1  Bacillus amyloliquefaciens                99.176277
2         Bacillus velezensis                99.176277
3         Bacillus velezensis                99.176277
4         Bacillus velezensis                99.176277


> Géneros a los que puede corresponder tu secuencia:


        - Bacillus


```


## Conclusión
Biopython es un set de herramientas que facilita el trabajo computacional en biología y el trabajo bioinformático. Con sus librerías especializadas en aplicaciones biológicas se puede eficientar las tareas bioinformáticas como buscar infomación en bases de datos, navegar un archivo o trabajar directamente con secuencias


## Referencias
Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. Journal of molecular biology, 215(3), 403-410.

Rodicio, M. del R., Mendoza, M, del C.(2004). Identificación bacteriana mediante secuenciación del ARNr 16S: fundamento, metodología y aplicaciones en microbiología clínica. Enfermedades Infecciosas y Microbiología Clínica, 22(4), 238.

Odom, A. R., Faits, T., Castro-Nallar, E., Crandall, K. A., & Johnson, W. E. (2023). Metagenomic profiling pipelines improve taxonomic classification for 16S amplicon sequencing data. Scientific reports, 13(1), 13957.
