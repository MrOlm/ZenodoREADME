# Data supporting publication: Ultra-deep Sequencing of Hadza Hunter-Gatherers Recovers Vanishing Gut Microbes

Link to Zenodo archive: https://doi.org/10.5281/zenodo.7782709

This upload contains several datasets to support the publication "Ultra-deep Sequencing of Hadza Hunter-Gatherers Recovers Vanishing Gut Microbes", including genomes, mapping databases, and other relevant data. The preprint publication is available on [bioRxiv](https://doi.org/10.1101/2022.03.30.486478). The citation for the peer-reviewed publication is not yet available as the time of this upload.

# Contents Directory

* **Assemblies_Proteins.zip** - A folder containing the called amino acid sequences for each metagenomic assembly, in .faa format. Generated using Prodigal in metagenome mode. The naming convention for each file is: sample name, two underscores (`__`), and then the name of the metagenomic assembly. You can find information linking each sample name to it's metadata are found in Supplemental Table S1 (within the file `Data_SupplementalTables.zip`).

* **Assemblies_ProteinAnnotation.zip** - A folder containing annotations for the called amino acid sequences for each metagenomic assembly. Annotations are against the `UniRef100` database and were generated using the program [tRep](https://github.com/MrOlm/tRep). Same naming convention as `Assemblies_Proteins.zip`.

* **Data_GeneDetection.zip** - A folder containing the gene-level inStrain profile data for all samples mapped to the HPlusDB_BacteriaArchaea database. The format of these files can be understood [at this link in the inStrain documentation](https://instrain.readthedocs.io/en/latest/example_output.html#gene-info-tsv). Same naming convention as `Assemblies_Proteins.zip`.

* **Data_GenomeDetection.zip** - A folder containing genome-level detection information. There are 3 files:
  * **BacArch_inStrainProfile_full.csv.gz** - full inStrain profile information for all samples mapped to the HPlusDB_BacteriaArchaea database. The format of these files can be understood [at this link in the inStrain documentation](https://instrain.readthedocs.io/en/latest/example_output.html#genome-info-tsv). Information on the genomes can be found in Supplemental Table S1 (within the file `Data_SupplementalTables.zip`).
  * **Bacteriphage_QuickProfile_breadth75.csv.gz** - inStrain quick-profile information for all samples mapped to the HPlusDB_Bacteriophage database. Information on the genomes can be found in Supplemental Table S1 (within the file `Data_SupplementalTables.zip`).
  * **Eukaryote_QuickProfile.csv.gz** - inStrain quick-profile information for all samples mapped to the HPlusDB_Eukaryotes database. Information on the genomes can be found in Supplemental Table S1 (within the file `Data_SupplementalTables.zip`).
  * **BachArch_GlobalMgx_QuickProfile.csv.gz** - inStrain quick-profile information for all public metagenome samples mapped to the BacArch database. Information on the genomes can be found in Supplemental Table S1 (within the file `Data_SupplementalTables.zip`) and information on the samples can be found in Supplementary Table S4.

* **Data_SupplementalTables.zip** - A folder containing copies of all Supplemental Tables associated with this manuscript. Descriptions are as follows:
  * **Supplementary Table 1** - Sample-level information for all metagenomes from the Hadza, Nepali, and Californian cohorts, along with links to raw reads and metagenomic assemblies.
  * **Supplementary Table 2** - Information about all genome recovered and analyzed in this study (including representative genomes and other genomes)
  * **Supplementary Table 3** - Roster of strains isolated from Hadza stool (including cultivation information)
  * **Supplementary Table 4** - Global metagenomics data set broken down by sample
  * **Supplementary Table 5** - Prevalence/abundance data for each species-level representative genome in our bacterial/archaeal species-level genome database
  * **Supplementary Table 6** - KO and Pfam gene level information
  * **Supplementary Table 7** - Strain sharing data between Hadza adult samples 

* **Genomes_BacteriaArchaea.zip** - A folder containing all Bacteria and Archaea genomes recovered in this study; 54,779 genomes in total. Supplemental Table S2 contains information on all recovered genomes.

* **Genomes_Bacteriophage.zip** - A folder containing all Bacteriophage genomes recovered in this study; 43,282 genomes in total. Supplemental Table S2 contains information on all recovered genomes 

* **Genomes_Eukaryotes.zip** - A folder with two subfolders: one containing all 21 Eukaryotic genomes recovered in this study, and another containing all protein .faa files for all Eukaryotic genomes. Supplemental Table S2 provides information on all recovered genomes

* **Genomes_Isolates.zip** - A folder containing all isolate genomes recovered in this study; 50 genomes in total. Supplemental Table S3 contains information on all isolate genomes.

* **HPlusDB_BacteriaArchaea.zip** - A folder containing the database of representative bacterial and archaeal genomes used in this study. Within this folder are the following files:
  * **HPlusDB_BacteriaArchaea.fasta.gz** - A .fasta file of all representitive genomes concatenated together. You can filter SupplementalTable_S2.csv to just retain rows where the column `Representitive_genome` is `TRUE` to get a list of genomes included in this fasta file.
  * **HPlusDB_BacteriaArchaea.stb.gz** - A [scaffold to bin file](https://instrain.readthedocs.io/en/latest/overview.html#term-scaffold-to-bin-file) denoting the bin to which each scaffold in HPlusDB_BacteriaArchaea.fasta belongs
  * **HPlusDB_BacteriaArchaea.genes.fna.gz** - A .fasta file of the nucleotide sequences of all genes called on all genomes in this database. Can be used as an input to inStrain.
  * **HPlusDB_BacteriaArchaea.genes.faa.gz** - A .faa file of a the animo acid sequences of all genes called on all genoems in this database. Can be used to annotate gennes against other databases
  * **bt2_index/** - A folder containing a bowtie2 index of HPlusDB_BacteriaArchaea.fasta. Used to map read with the program Bowtie2
  * **annotation/** - A folder containing annotations for HPlusDB_BacteriaArchaea.genes.faa against 5 different databases- Pfam, UniRef100, KOfamscan, abricate, and dbCAN

* **HPlusDB_Bacteriophage.zip** - A folder containing the database of representative bacteriophage genomes used in this study. Within this folder are the following files:
  * **HPlus_BacteriophageDB.fasta.gz** - A .fasta file of all representitive genomes concatenated together
  * **HPlus_BacteriophageDB.stb.gz** - A [scaffold to bin file](https://instrain.readthedocs.io/en/latest/overview.html#term-scaffold-to-bin-file) denoting the bin to which each scaffold in HPlus_BacteriophageDB.fasta belongs
  * * **bt2_index/** - A folder containing a bowtie2 index of HPlus_BacteriophageDB.fasta. Used to map reads with the program Bowtie2

* **HPlusDB_Eukaryotes.zip** - A folder containing the database of representative eukaryotic genomes used in this study. Within this folder are the following files:
  * **HPlusDB_Eukaryotes.fasta** - A .fasta file of all representitive genomes concatenated together
  * **HPlusDB_Eukaryotes.stb** - A [scaffold to bin file](https://instrain.readthedocs.io/en/latest/overview.html#term-scaffold-to-bin-file) denoting the bin to which each scaffold in HPlusDB_Eukaryotes.fasta belongs
  * * **bt2_index/** - A folder containing a bowtie2 index of HPlusDB_Eukaryotes.fasta. Used to map reads with the program Bowtie2

# Downloading additional raw data

In addition to the data in the repository, all raw reads and assemblies are hosted on NCBI's servers. Links to these files are available in Supplemental_Table_S1.csv, within the file `Data_SupplementalTables.zip`. Below are examples of how to easily download these files.

## Download raw reads

The column `Reads_accession` in the file `Supplemental_Table_S1.csv` can be used to download the metagenomic reads for all metagenomes produced in this study. This is best done using the [sra-tools package](https://github.com/ncbi/sra-tools), which can be installed with conda using the command `conda install -c bioconda sra-tools`.

Reads for individual samples can be downloaded using the sample's `Reads_accession` with a command like the following:

    fastq-dump --split-files ERR7738514

The following bash loop can also be used to download reads for all the metagenomes produced in this study:

    cat SupplementalTable_S1.csv | tail -n +2 | while read line; do   SRA=$(echo "$line" | cut -d ',' -f 22); echo  fastq-dump --split-files "$SRA"; done

## Download raw assemblies

Ftp links to download the raw assemblies are also available in the file `Supplemental_Table_S1.csv` in the column `Assembly_link`. These can be easily downloaded using a command like the following:

    wget ftp://ftp.sra.ebi.ac.uk/vol1/sequence/ERZ456/ERZ4567069/contig.fa.gz -O Pilot_MoBio_Fiber-Hadza-Nepal_C_23_7024.assembly.fasta.gz
    
The following bash loop can also be used to download reads for all the metagenomic assemblies produced in this study:

    cat SupplementalTable_S1.csv | tail -n +2 | while read line; do sample_name=$(echo "$line" | awk -F ',' '{print $1}'); assembly_link=$(echo "$line" | awk -F ',' '{print $24}'); wget "$assembly_link" -O "${sample_name}.assembly.fasta.gz"; done

# Tutorials

Below are some breif tutorials on how to use the above data to perform common tasks.

## Analyze your reads with the genomes generated in this study

This tutorial uses Bowtie2 and inStrain to analyze YOUR reads against the genoems in THIS database. It is similar to [this online tutorial](https://instrain.readthedocs.io/en/latest/tutorial.html#tutorial-2-running-instrain-using-a-public-genome-database)

If you don't have reads but would like to follow along, you can download the smallest metagenome generated in this study using the following command:

    fastq-dump --split-files ERR7737607
    
The first step is to map reads to the genome database. The Bowtie2 index is available for download in the file  `HPlusDB_BacteriaArchaea.zip`. Once downloaded you need to decompress all files in the bt2_index folder using a program like `gzip`. You can then map the reads to genearte a .sam file using a command like the following:


    bowtie2 -p 10 -x HPlusDB_BacteriaArchaea.bt2 -1 ERR7737607_1.fastq -2 ERR7737607_2.fastq > HPlusDB_BacteriaArchaea-vs-ERR7737607.sam

Once the mapping is done you can run inStrain to analyze the resulting .sam file using verious other files profiled in the `HPlusDB_BacteriaArchaea.zip` file. A example full inStrain profile command would look like the following:

    inStrain profile HPlusDB_BacteriaArchaea-vs-ERR7737607.sam HPlusDB_BacteriaArchaea.fasta -o HPlusDB_BacteriaArchaea-vs-ERR7737607.IS -p 10 -g HPlusDB_BacteriaArchaea.genes.fna -s HPlusDB_BacteriaArchaea.stb --database_mode
    
See the [inStrain documentation](https://instrain.readthedocs.io/en/latest/example_output.html#instrain-profile) for descriptions of the resultant output files.

The file `SupplementalTable_S2.csv` has the taxonomy for all genomes in this database. You can filter this table to just retain rows where the column `Representitive_genome` is `TRUE` to get a list of genomes included in the HPlusDB_BacteriaArchaea database. The folder `annotation` in the file `HPlusDB_BacteriaArchaea.zip` has annotations for all profiled genes.

[See full inStrain documentation for more information on this technique](https://instrain.readthedocs.io/en/latest/index.html)

## Integrate the genomes in your study with the genomes generated in this study

The above tutorial has you mapping your reads to the primary genome database using in this study. This tutorial shows how to make a **new** mapping database by merging a new set of genomes with genome database used in this study. It's based on [this inStrain tutorial](https://instrain.readthedocs.io/en/latest/tutorial.html#tutorial-3-merging-custom-genomes-with-an-existing-genome-database).

The first step is to create individual files for each genome in the `HPlusDB_BacteriaArchaea` database. We can do this in two ways.

### Method 1) Create individual genome files using the .stb file

This method uses the [parse_stb.py script](https://github.com/MrOlm/drep/blob/master/helper_scripts/parse_stb.py) that comes with the program [dRep](https://drep.readthedocs.io/en/latest/#). Running the follings commands:

  mkdir HPlus_BacteriaArchaea_genomes/
  
  parse_stb.py -f HPlusDB_BacteriaArchaea.fasta.gz	 -s HPlusDB_BacteriaArchaea.stb -o HPlus_BacteriaArchaea_genomes/
  
will generate a individual .fasta files for all genomes in the folder HPlus_BacteriaArchaea_genomes/

### Method 2) Copy individual genomes from `Genomes_BacteriaArchaea.zip` and UHGG

Within the `Genomes_BacteriaArchaea.zip` folder are all genomes recovered in this study. Filtering **Supplemental_Table_S2.excel** to just the rows where `Representitive_genome` is `TRUE` will give you a list of the genomes that are in the `HPlusDB_BacteriaArchaea` database. You can then copy over all of the genomes recovered from `Genomes_BacteriaArchaea.zip` to a new folder, and then copy all of representitive genomes from `UHGG` into that same folder.

Both **Method 1** and **Method 2** will result in you having a folder of the individual genomes making up `HPlusDB_BacteriaArchaea`. To dereplicate these genomes with your genomes, we can again use the program dRep with a command like the following:

  dRep dereplicate MergedGenomeSet -g HPlus_BacteriaArchaea_genomes/* custom_genome/* –S_algorithm fastANI –multiround_primary_clustering –clusterAlg greedy -ms 10000 -pa 0.9 -sa 0.95 -nc 0.30 -cm larger -p 16
  
The specific thresholds you use for dereplication are important and some thoughts about choosing thresholds is available at [this link](https://instrain.readthedocs.io/en/latest/important_concepts.html).

Once this dRep run is finished, you can follow the instructions at the top of [this tutorial](https://instrain.readthedocs.io/en/latest/tutorial.html#tutorial-2-running-instrain-using-a-public-genome-database) to create a new mapping database for this merged genome database. To prioritize your custom genomes over the database genomes, you can use the flag `extra_weight_table` within dRep.

[See full dRep documentation for more information on this technique](https://drep.readthedocs.io/en/latest/)
