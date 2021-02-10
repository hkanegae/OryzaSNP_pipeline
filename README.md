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
The Plant Genome. 2020;e20005. 
[doi.org/10.1002/tpg2.20005](https://doi.org/10.1002/tpg2.20005)  
[Pipeline LINK](https://github.com/hkanegae/OryzaSNP_pipeline/blob/master/PMID30458025.md) 

### Coupling day length data and genomic prediction tools for predicting time-related traits under complex scenarios
Diego Jarquin, Hiromi Kajiya-Kanegae, Chen Taishen, Shiori Yabe, Reyna Persa, Jianming Yu, Hiroshi Nakagawa, Masanori Yamasaki, Hiroyoshi Iwata
Sci Rep. 2020 A10:13382. 
[doi: 10.1038/s41598-020-70267-9.](https://www.nature.com/articles/s41598-020-70267-9)  
PMID: 32770083   
[Pipeline LINK](https://github.com/hkanegae/OryzaSNP_pipeline/blob/master/PMID30458025.md) 
***
### Predicting Rice Heading Date Using an Integrated Approach Combining a Machine Learning Method and a Crop Growth Model
Tai-Shen Chen, Toru Aoike, Masanori Yamasaki, Hiromi Kajiya-Kanegae and Hiroyoshi Iwata 
[doi: 10.3389/fgene.2020.599510](https://www.frontiersin.org/articles/10.3389/fgene.2020.599510/full)  
PMID: 33391352   
[Pipeline LINK](https://github.com/hkanegae/OryzaSNP_pipeline/blob/master/PMID30458025.md) 
***

## Analysis workflow for detection of genome-wide variations in TASUKE+ of RAP-DB was modified.

### REF=RGSP-1.0 genome (including organella and unanchored contig sequences) [FASTA](https://rapdb.dna.affrc.go.jp/download/archive/genome-wide_variations/IRGSP-1.0_genome_M_C_unanchored.fa.gz)

### trimmomatic (v0.38)
$ java -jar trimmomatic-0.38.jar PE \
    -phred33 read.r1.fastq.gz read.r2.fastq.gz \
    read.pe.r1.fastq.gz read.se.r1.fastq.gz read.pe.r2.fastq.gz read.se.r2.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:30

### Mapping / bwa v0.7.17
$ java -jar picard.jar FastqToSam \
    FASTQ=read.pe.r1.fastq.gz \
    FASTQ2=read.pe.r2.fastq.gz \
    OUTPUT=uBAM.bam \
    READ_GROUP_NAME=${SAMPLE_ID} \
    SAMPLE_NAME=${SAMPLE_ID} \
    LIBRARY_NAME=${SAMPLE_ID} \

$ java -jar picard.jar MergeBamAlignment \
    ALIGNED=alignment.sort.bam \
    UNMAPPED=uBAM.bam \
    OUTPUT=alignment.merge.bam \
    REFERENCE_SEQUENCE=genome.fa

### samtools flagstat / samtools (v1.9)
$samtools flagstat "$dir2"/"$name".sorted.bam > "$dir3"/"$name"_flagstat.txt

### FixMate information / picard-tools-2.5.0
java -jar picard.jar FixMateInformation I="$dir2"/"$name".sorted.bam O="$dir5"/"$name".fxmt.bam SO=coordinate CREATE_INDEX=TRUE

### Remove PCR duplicates  / Picard (v2.18.17)
$ java -jar picard.jar MarkDuplicates \
    INPUT=alignment.merge.bam \
    OUTPUT=alignment.rmdup.bam \
    METRICS_FILE=rmdup.matrix \
    REMOVE_DUPLICATES=true \
    MAX_RECORDS_IN_RAM=1000000 \
    TMP_DIR=./tmp
$ samtools index alignment.rmdup.bam

### Variant detection and filtering by GATK / GATK (v4.0.11.0)
$ gatk HaplotypeCaller \
    --input alignment.rmdup.bam \
    --output variants.g.vcf.gz \
    --reference genome.fa \
    -max-alternate-alleles 2 \
    --emit-ref-confidence GVCF

### GenomicsDBImport
array=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12)

for name in ${array[@]}; do echo $name; 
gatk --java-options "$JAVA_MEM" GenomicsDBImport --reference genome.fa -V variants_A.g.vcf.gz -V variants_B.g.vcf.gz -V variants_C.g.vcf.gz --genomicsdb-workspace-path variants_"$name" --intervals "$name"
done

 
### GenotypeGVCFs
for name in ${array[@]}; do echo $name;
gatk --java-options "$JAVA_MEM" GenotypeGVCFs --reference genome.fa -V gendb://variants_"$name" -G StandardAnnotation --new-qual -O variants.genotype_"$name".vcf.gz
done

### VariantFiltration

for name in ${array[@]}; do echo $name;
gatk --java-options "$JAVA_MEM" VariantFiltration --reference genome.fa --variant variants.genotype_"$name".vcf.gz --output variants.filter.genotype_"$name".vcf.gz --filter-expression "QD < 5.0 || FS > 50.0 || SOR > 3.0 || MQ < 50.0 || MQRankSum < -2.5 || ReadPosRankSum < -1.0 || ReadPosRankSum > 3.5" --filter-name "FILTER"
done

### SelectVariants

for name in ${array[@]}; do echo $name;
gatk --java-options "$JAVA_MEM" SelectVariants --reference genome.fa --variant variants.filter.genotype_"$name".vcf.gz  --output variants.varonly.vcf.gz  --exclude-filtered --select-type-to-include SNP  --select-type-to-include INDEL
done

### merge varonly.vcf
java $JAVA_MEM -jar $PICARD_HOME/picard.jar MergeVcfs O=variants.varonly.vcf.gz I= variants_chr01.varonly.vcf.gz I=variants_chr02.varonly.vcf.gz I= variants_chr03.varonly.vcf.gz I=variants_chr04.varonly.vcf.gz I= variants_chr05.varonly.vcf.gz I=variants_chr06.varonly.vcf.gz I= variants_chr07.varonly.vcf.gz I=variants_chr08.varonly.vcf.gz I= variants_chr09.varonly.vcf.gz I=variants_chr10.varonly.vcf.gz I= variants_chr11.varonly.vcf.gz I= variants_chr12.varonly.vcf.gz


