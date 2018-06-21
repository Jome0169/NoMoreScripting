## Annotation Scripts for Naming, and Updating a Genome Annotation

Genomic annotations can be increidble tedious and challenging projects to
undergo. With poor conventions, and nearly no universal standards, it can be
hard to identify what makes a "well done" annotation, in format. The scripts in
this directoty attempt to address this problem for the Kingdom of Plantae.
Final functional annotations created by this series of scripts abide by the
naming convention and standards of Arabidopsis Thaliana laid out by the TAIR
consortium.

To function correctly this program expects a MAKER output. Also -- NOTE: Since
gene naming is done by the coordinates of a gene on a chromsome, it behooves
you to ensure you cannot add any additional genes to the annotation before
embarking on this endheavor. As adding new genes will either need to be done by
hand, or require rerunning the entire pipeline over again.

A quick guide to some basic precursor steps can be found over at <Link>


### Step 1: Create a mapping file

In order to create a mapping file to change gene names - you must first sort
the gff3 output files using gfftools. If you do not sort, gene names will be
out of place and incorrect.

Before creating a mapping file, ensure that you have an adequete descriptive
prefix that will be you geneoms unique GeneID prefix. Generally its usually the
first two letters of the genus - followed by the first two letters of the
species.

```
python3 mapping_file.py -g Final.total.sorted.gff3 -pre Cuco -o mapping.file.txt

```

Your mapping file should look something like:
```
maker-Chr1-snap-gene-0.223	Cucu1G000010
maker-Chr1-snap-gene-0.223-mRNA-1	Cucu1G000010.1
maker-Chr1-augustus-gene-0.133	Cucu1G000020
maker-Chr1-augustus-gene-0.133-mRNA-1	Cucu1G000020.1
maker-Chr1-snap-gene-0.226	Cucu1G000030
maker-Chr1-snap-gene-0.226-mRNA-1	Cucu1G000030.1
maker-Chr1-augustus-gene-0.135	Cucu1G000040
maker-Chr1-augustus-gene-0.135-mRNA-1	Cucu1G000040.1
augustus_masked-Chr1-processed-gene-0.54	Cucu1G000050
augustus_masked-Chr1-processed-gene-0.54-mRNA-1	Cucu1G000050.1
augustus_masked-Chr1-processed-gene-0.55	Cucu1G000060
augustus_masked-Chr1-processed-gene-0.55-mRNA-1	Cucu1G000060.1
maker-Chr1-augustus-gene-0.136	Cucu1G000070
maker-Chr1-augustus-gene-0.136-mRNA-1	Cucu1G000070.1
maker-Chr1-snap-gene-0.195	Cucu1G000080
maker-Chr1-snap-gene-0.195-mRNA-1	Cucu1G000080.1
maker-Chr1-augustus-gene-0.174	Cucu1G000090
maker-Chr1-augustus-gene-0.174-mRNA-1	Cucu1G000090.1
maker-Chr1-snap-gene-0.197	Cucu1G000100
maker-Chr1-snap-gene-0.197-mRNA-1	Cucu1G000100.1
maker-Chr1-augustus-gene-0.142	Cucu1G000110
```

NOTE: If you do not have an equal number of genes and mRNAs in your gff3 file,
this script will give you an INCORRECT mapping file. So before you go any
further double check your original gff3 file with:

```
awk '{print $3}' Final.genes.gff | grep 'mRNA' | wc
awk '{print $3}' Final.genes.gff | grep 'gene' | wc
```

Both numbers output should be identical. Otherwise you may need to go through
and remove alternative splicing events that made their way into your finalized
maker gff3.


### Step 2: Rename GFF3 File With New Names

Now that we have a mapping file it's a rather quick fix to rename all genes in
both the Fasta and gff3 file. To rename the gff and fasta files, the commands
are: 

```
python3 ~/Programming/NoMoreScripting/nomorescripting/fasta/rename_fasta_file.py -f Final.protein.rawnames.fasta -m mapping.file.txt > corrected.protein.fast

python3 ~/Programming/NoMoreScripting/nomorescripting/annotation/rename_gff_file.py -g Final.genes.gff -m mapping.file.txt > corrected.gff
```

### Step 3: Run Functional Annotation Software

Now, in order to actually identify what genes actually have homology with other
genes, it's essential to blast the renamed protein file to multiple databases
of similar species, as well as run various softwares such as blast2go, AHRD,
and interproscan.

Once these softwares have all been run and their various output files
collected, it's then possible to add the functional annotation component to
each gene in the Fasta files, as well as attach their functional information to
the gff3 file as well.

### Step 4: Adding Functional Information







