# Folder and data preparation
## Create directory structure for annotation, genome data, Salmon index, sample processing, trimming
mkdir -p /data/{annotation,genome,index,samples,trim}
mkdir -p /data/results/rMATS/res/tmp
mkdir -p /data/results/vast_tools/{compare/tidy,compare_expr,diff/tidy}

## Download the GTF annotation file for Mus musculus (version: 106) directly to /data/annotation directory (ENSEMBL) and unzip the downloaded GTF file
wget -P /data/annotation "http://ftp.ensembl.org/pub/release-106/gtf/mus_musculus/Mus_musculus.GRCm39.106.chr.gtf.gz" && gunzip /data/annotation/Mus_musculus.GRCm39.106.chr.gtf.gz

## Download the primary assembly, cDNA and ncRNA FASTA files for Mus musculus (version: 106) ENSEMBL repository
wget -P /data/genome "http://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
wget -P /data/genome "http://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz"
wget -P /data/genome "http://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" && gunzip /data/genome/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

## Concatenate cDNA and ncRNA FASTA files into a single file
zcat /data/genome/Mus_musculus.GRCm39.cdna.all.fa.gz /data/genome/Mus_musculus.GRCm39.ncrna.fa.gz > /data/genome/Mus_musculus.GRCm39.cdna.ncrna.fa

## Transfer raw FASTQ files from local machine to server (anonymized)
scp /path/to/local/files/* user@ip_server:/data/samples/

## Combine all the raw sequencing files lines into one file for different replicates
cat /data/samples/condition?_loc_R1_L00?* > /data/samples/condition?_R1.fastq.gz
cat /data/samples/condition?_loc_R2_L00?* > /data/samples/condition?_R2.fastq.gz

cat /data/samples/wt?_loc_R1_L00?* > /data/samples/wt?_R1.fastq.gz
cat /data/samples/wt?_loc_R2_L00?* > /data/samples/wt?_R2.fastq.gz

## Remove the intermediate files (raw names) to keep the directory clean
rm /data/samples/condition?_loc_R?_L00?* 
rm /data/samples/wt?_loc_R?_L00?*

## Loop through paired-end FASTQ files and process them with trim_galore (version: 0.6.4_dev)
for i in $(ls /data/samples/*.fastq.gz | xargs -n 2 basename | cut -d '_' -f 1); do trim_galore --paired -j 7 --dont_gzip /data/samples/${i}_R1.fastq.gz /data/samples/${i}_R2.fastq.gz -o /data/results/ --fastqc; done &

# rMATS (version: v4.1.2)
## Software download in /data/rMats/rmats_turbo_v4_1_2/. We add the PATH into /etc/enviroment to be available for all users using the server.
## We create a metadata file /data/samples/txt/wt.txt and /data/samples/txt/condition.txt
nano /data/samples/txt/wt.txt
cat /data/samples/txt/condition.txt

## STAR mapping (sjdbOverhang 149, because read length is 150 bp)
STAR --runMode genomeGenerate --genomeDir /data/index/STAR/genome_mm --genomeFastaFiles /data/genome/Mus_musculus.GRCm39.dna.primary_assembly.fa --sjdbGTFfile  /data/annotation/Mus_musculus.GRCm39.106.chr.gtf --runThreadN 50 --sjdbOverhang 149

## to open many files at once
ulimit -n 3076

## rMATS works with Python 2.7 instead of 3, we create a conda enviroment for that
conda activate python2.7

## Run the rMATS differential splicing analysis using Python
## - --s1 specifies the input file for wild-type samples (wt.txt)
## - --s2 specifies the input file for condition samples (condition.txt)
## - --gtf provides the GTF annotation file for the Mus musculus genome (GRCm39)
## - --bi specifies the path to the STAR genome index
## - -t paired indicates paired-end RNA-Seq data
## - --readLength 150 sets the read length for the RNA-Seq data to 150 base pairs
## - --cstat 0.05 sets a 5% cutoff for statistical significance
## - --nthread 30 uses 30 CPU threads for parallel processing
## - --od specifies the directory for the output results
## - --tmp specifies the directory for temporary files
## - & runs the command in the background to allow other tasks to continue
python rmats.py --s1 /data/samples/txt/wt.txt --s2 /data/samples/txt/condition.txt --gtf /data/annotation/Mus_musculus.GRCm39.106.chr.gtf --bi /data/index/STAR/genome_mm/ -t paired --readLength 150 --cstat 0.05 --nthread 30 --od /data/results/rMATS/res/ --tmp /data/results/rMATS/res/tmp/ &

# Sashimiplot (version: 2.0.4)
## This section generates sashimi plots for different types of alternative splicing events.
for i in SE MXE A5SS A3SS RI; do rmats2sashimiplot --b1 ./wt.bam --b2 ./condition.bam --l1 wt --l2 condition -t ${i} -e /data/results/rMATS/res/${i}.MATS.JC.txt -o /data/results/rMATS/res/rmatstosashimiplot/JC/${i}; done &

# vast-tools (version: v2.5.1)
## Align reads using vast-tools for each sample.
for i in $(ls /data/samples/trim/*.fq | xargs -n 2 basename | cut -d '_' -f 1); do /bin/vast-tools/vast-tools align -c 30 /data/samples/trim/${i}_val_1.fq /data/samples/trim/${i}_val_2.fq -sp mm10 -n ${i} -o /data/results/vast_tools/ --expr; done &

## Combine the alignment results for all samples.
/bin/vast-tools/vast-tools combine -c 6 -o /data/results/vast_tools/ -C -sp mm10

## Compare inclusion levels between two conditions.
/bin/vast-tools/vast-tools compare /data/results/vast_tools/INCLUSION_LEVELS_FULL-mm10-2.tab -a condition -b wt --name_A condition --name_B wt > /data/results/vast_tools/compare/compare.txt

## compare > tidy
## Tidy up the comparison output for better readability.
/bin/vast-tools/vast-tools tidy /data/results/vast_tools/DiffAS-mm10-2-dPSI15-range5-min_ALT_use25-upreg_ALT_condition-vs-wt.tab --samples condition,wt --add_names --log > /data/results/vast_tools/compare/tidy/compare_tidy.txt

## Generate plots for the differential alternative splicing results.
/bin/vast-tools/vast-tools plot /data/results/vast_tools/DiffAS-mm10-2-dPSI15-range5-min_ALT_use25-upreg_ALT_condition-vs-wt.tab

## Compare gene expression levels between two conditions.
/bin/vast-tools/vast-tools compare_expr /data/results/vast_tools/cRPKM_AND_COUNTS-mm10-2.tab -a condition -b wt --name_A condition --name_B wt > /data/results/vast_tools/compare_expr/compare_expr.txt

## Perform differential splicing analysis at different thresholds.
/bin/vast-tools/vast-tools diff -a condition -b wt --sampleNameA=Condition --sampleNameB=WT -o /data/results/vast_tools/diff/
/bin/vast-tools/vast-tools diff -m 0.2 -d m02 -a condition -b wt --sampleNameA=Condition --sampleNameB=WT -o /data/results/vast_tools/diff/
/bin/vast-tools/vast-tools diff -m 0.4 -d m04 -a condition -b wt --sampleNameA=Condition --sampleNameB=WT -o /data/results/vast_tools/diff/

## Tidy up the differential splicing output.
/bin/vast-tools/vast-tools tidy INCLUSION_LEVELS_FULL-mm10-2.DIFF.txt --samples condition,wt --add_names --log > /data/results/vast_tools/diff/tidy/output_diff_tidy.txt

## Tidy up the inclusion levels output.
/bin/vast-tools/vast-tools tidy INCLUSION_LEVELS_FULL-mm10-2.tab --samples condition,wt --add_names --log > /data/results/vast_tools/diff/tidy/output_tidy.txt

