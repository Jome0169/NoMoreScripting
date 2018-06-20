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


### Step 2: Rename GF















### 


