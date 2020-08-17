## OryzaSNP_pipeline
A description of the analysis step of Oryza sativa variant calling pipeline.

### Description of grain weight distribution leading to genomic selection for grain-filling characteristics in rice.  
Yabe S, Yoshida H, Kajiya-Kanegae H, Yamasaki M, Iwata H, Ebana K, Hayashi T, Nakagawa H.   
PLoS One. 2018 Nov 20;13(11):e0207627.  
[doi: 10.1371/journal.pone.0207627.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0207627)
PMID: 30458025   
[Pipeline LINK](https://github.com/hkanegae/OryzaSNP_pipeline/blob/master/PMID30458025.md)  

### Choosing the optimal population for a genome‐wide association study: A simulation of whole‐genome sequences from rice
Kosuke Hamazaki  Hiromi Kajiya‐Kanegae  Masanori Yamasaki  Kaworu Ebana  Shiori Yabe  Hiroshi Nakagawa  Hiroyoshi Iwata
The Plant Genome. 2019;e20005. 
[doi.org/10.1002/tpg2.20005](https://doi.org/10.1002/tpg2.20005)

### Coupling day length data and genomic prediction tools for predicting time-related traits under complex scenarios
Diego Jarquin, Hiromi Kajiya-Kanegae, Chen Taishen, Shiori Yabe, Reyna Persa, Jianming Yu, Hiroshi Nakagawa, Masanori Yamasaki, Hiroyoshi Iwata
Sci Rep. 2020 A10:13382. 
[doi: 10.1038/s41598-020-70267-9.](https://www.nature.com/articles/s41598-020-70267-9)
***

### REF=IRGSP-1.0_genome.fasta
### trimmomatic
java -jar trimmomatic-0.36.jar PE -threads 4 -phred33 "$dir"/"$name"_1.fastq.gz "$dir"/"$name"_2.fastq.gz "$dir"/"$name"_1_paired.fastq.gz "$dir"/"$name"_1_unpaired.fastq.gz "$dir"/"$name"_2_paired.fastq.gz "$dir"/"$name"_2_unpaired.fastq.gz ILLUMINACLIP:/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

### Mapping / bwa-0.7.12
$bwa mem -M -t 4 oryza_index "$dir"/"$name"_1_paired.fastq.gz "$dir"/"$name"_2_paired.fastq.gz | $samtools view -bS - | $samtools sort -T tmpsam"$name" -@4 -o "$dir2"/"$name".sorted.bam

### samtools flagstat / samtools-1.3.1
$samtools flagstat "$dir2"/"$name".sorted.bam > "$dir3"/"$name"_flagstat.txt

### sort_bam index
$samtools index "$dir2"/"$name".sorted.bam

### FixMate information / picard-tools-2.5.0
java -jar picard.jar FixMateInformation I="$dir2"/"$name".sorted.bam O="$dir5"/"$name".fxmt.bam SO=coordinate CREATE_INDEX=TRUE

### Mark duplicate reads
java -jar picard.jar MarkDuplicates I="$dir5"/"$name".fxmt.bam O="$dir5"/"$name".mkdup.bam M="$dir5"/"$name".metrics.txt CREATE_INDEX=TRUE

### Add or replace read groups
java -jar picard.jar AddOrReplaceReadGroups I="$dir5"/"$name".mkdup.bam O="$dir5"/"$name".addrep.bam RGPL=illumina RGLB=lib1 RGPU=unit1 RGSM="$name" RGID="$name"

$samtools index "$dir5"/"$name".addrep.bam

### GATK RealignerTargetCreator / GATK 3.7-0-gcfedb67
java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF -I "$dir5"/"$name".addrep.bam -o "$dir5"/"$name".addrep.intervals

### GATK IndelRealigner
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R $REF -I "$dir5"/"$name".addrep.bam -targetIntervals "$dir5"/"$name".addrep.intervals -o "$dir5"/"$name".realn.bam

### CollectRawWgsMetrics
java -jar picard.jar CollectRawWgsMetrics I="$dir5"/"$name".realn.bam O="$dir6"/"$name"_metrics.txt R=$REF INCLUDE_BQ_HISTOGRAM=true

### SNP call

java -Xmx8g -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -I "$dir5"/"$name".realn.bam -R $REF -glm BOTH -gt_mode DISCOVERY -o $name.vcf -nt 4
