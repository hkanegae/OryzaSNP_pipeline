## OryzaSNP_pipeline
A description of the analysis step of Oryza sativa variant calling pipeline.

### REF=IRGSP-1.0_genome.fasta
### trimmomatic
java -jar trimmomatic-0.36.jar PE -threads 4 -phred33 "$dir"/"$name"_1.fastq.gz "$dir"/"$name"_2.fastq.gz "$dir"/"$name"_1_paired.fastq.gz "$dir"/"$name"_1_unpaired.fastq.gz "$dir"/"$name"_2_paired.fastq.gz "$dir"/"$name"_2_unpaired.fastq.gz ILLUMINACLIP:/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

### Mapping / bwa-0.7.12
$bwa mem -M -t 4 oryza_index "$dir"/"$name"_1_paired.fastq.gz "$dir"/"$name"_2_paired.fastq.gz | $samtools view -bS - | $samtools sort -T tmpsam"$name" -@4 -o "$dir2"/"$name".sorted.bam

### samtools fragstat / samtools-1.3.1
$samtools flagstat "$dir2"/"$name".sorted.bam > "$dir3"/"$name"_flagstat.txt

### sort_bam index
$samtools index "$dir2"/"$name".sorted.bam

### FixMate information / picard-tools-2.5.0
java -jar picard.jar FixMateInformation I="$dir2"/"$name".sorted.bam O="$dir5"/"$name".fxmt.bam SO=coordinate CREATE_INDEX=TRUE

### Mark duplicate reads
java -jar picard.jar MarkDuplicates I="$dir5"/"$name".fxmt.bam O="$dir5"/"$name".mkdup.bam M="$dir5"/"$name".metrics.txt CREATE_INDEX=TRUE

### Add or replace read groups
java -jar picard.jar AddOrReplaceReadGroups I="$dir5"/"$name".mkdup.bam O="$dir5"/"$name".addrep.bam RGPL=illumina RGLB=lib1 RGPU=unit1 RGSM=20
$samtools index "$dir5"/"$name".addrep.bam

### GATK RealignerTargetCreator
java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF -I "$dir5"/"$name".addrep.bam -o "$dir5"/"$name".addrep.intervals

### GATK IndelRealigner
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R $REF -I "$dir5"/"$name".addrep.bam -targetIntervals "$dir5"/"$name".addrep.intervals -o "$dir5"/"$name".realn.bam

### CollectRawWgsMetrics
java -jar picard.jar CollectRawWgsMetrics I="$dir5"/"$name".realn.bam O="$dir6"/"$name"_metrics.txt R=$REF INCLUDE_BQ_HISTOGRAM=true

### SNP call

java -Xmx8g -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -I "$dir5"/"$name".realn.bam -R $REF -glm BOTH -gt_mode DISCOVERY -o $name.vcf -nt 4
