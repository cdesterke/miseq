#!/bin/bash
##miseq pipeline
##Christophe Desterke_nov2018
## sh miseq.sh R1.fastq R2.fastq 

nom_fichier=$(echo $1 | sed -re 's/(.*).bam/\1/')

## aligned fastq with bowti2 and gap parameters
bowtie2 -x genome -1 $1 -2 $2 -S aligned.sam --rdg 5,1 --rfg 5,1 --dpad 100
echo "--------------------"
echo "step1 : bowtie2 alignement"
UPTIME1=`uptime`

## bam compression and sorted with samtools
samtools view -bSu aligned.sam  | samtools sort - aligned.sorted
echo "--------------------"
echo "step2 : bam compression and sorted samtools"
UPTIME2=`uptime`

##picardtools remove duplicates
picard-tools MarkDuplicates INPUT=aligned.sorted.bam  OUTPUT=marked.bam METRICS_FILE=metrics CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true
echo "--------------------"
echo "step3 : picard remove duplicates"
UPTIME3=`uptime`

## picard ajouter un ReadGroup pour rendre compatible GATK
picard-tools AddOrReplaceReadGroups I=marked.bam O=markedRG.bam LB=whatever PL=illumina PU=whatever SM=whatever
echo "--------------------"
echo "step4 : picard add read group"
UPTIME4=`uptime`

## restriction of the bam to the myeloid panel
bedtools intersect -abam markedRG.bam -b sophia.bed > myeloid.bam
picard-tools ReorderSam I=markedRG.bam O=markedRGorder.bam REFERENCE=genome.fa
echo "--------------------"
echo "step5 : bedtools myeloid restriction"
UPTIME5=`uptime`

## samtools filtration Q30 quality
samtools view -h -b -q 30 myeloid.bam -o myeloidQ30.bam 
echo "--------------------"
echo "step6 : quality Q30 filtration"
UPTIME6=`uptime`

## samtools index of the bam file
samtools index myeloidQ30.bam
echo "--------------------"
echo "step7 : index of the bam"
UPTIME7=`uptime`


rm aligned.sam
rm aligned.sorted.bam
rm markedRG.bam
rm marked.bam
rm myeloid.bam


## GATK "local realignment around indels"

# creation d'une table des possibles indels
java -Xmx4g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -known Mills_and_1000G_gold_standard.indels.hg19.vcf -R genomeHG19.fa -o marked.bam.list -I myeloidQ30.bam --filter_reads_with_N_cigar

# realigne les reads autour de ces possibles indels
java -Xmx4g -Djava.io.tmpdir=/tmp -jar GenomeAnalysisTK.jar -I myeloidQ30.bam -R genomeHG19.fa -T IndelRealigner -known Mills_and_1000G_gold_standard.indels.hg19.vcf -targetIntervals marked.bam.list -o marked.realigned.bam --filter_reads_with_N_cigar

echo "--------------------"
echo "step8 : GATK realignement around indels"
UPTIME8=`uptime`

## picard fixmate (fixation de l'information mate comme l'alignement peut etre modifié pendant le processus de realignement)
picard-tools FixMateInformation I=marked.realigned.bam O=fixed.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

echo "--------------------"
echo "step9 : GATK fixmate"
UPTIME9=`uptime`

## base recalibrator
java -Xmx4G -jar GenomeAnalysisTK.jar -T BaseRecalibrator -I fixed.bam -R genomeHG19.fa -knownSites dbsnp.vcf -o recaldata.table

java -Xmx4G -jar GenomeAnalysisTK.jar -T PrintReads -R genomeHG19.fa -I fixed.bam -BQSR recaldata.table -o recal.bam

echo "--------------------"
echo "step10 : GATK base recalibration"
UPTIME10=`uptime`

## samtools 
samtools sort recal.bam recal.srt
samtools index recal.srt.bam

echo "--------------------"
echo "step11 : Samtools sort and index"
UPTIME11=`uptime`

## GATK haplotype caller 
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R genomeHG19.fa -I recal.srt.bam -o output.vcf 

echo "--------------------"
echo "step12 : HC calling"
UPTIME12=`uptime`

mv recal.srt.bam $(echo recal.srt.bam | sed "s/\./".$nom_fichier"\./")
mv recal.srt.bam.bai $(echo recal.srt.bam.bai | sed "s/\./".$nom_fichier"\./")




rm myeloidQ30.bam
rm myeloidQ30.bam.bai
rm marked.bam.list
rm marked.realigned.bam
rm fixed.bam
rm recaldata.table
rm recal.bam
rm recal.bai



##conversion du VCF en AVinput ANNOVAR
perl convert2annovar.pl -format vcf4  output.vcf -outfile outav

UPTIME13=`uptime`

perl table_annovar.pl outav humandb/ -buildver hg19 -out refgene -remove -protocol refGene,ensGene,cytoBand,genomicSuperDups,popfreq_all,esp6500si_all,1000g2012apr_all,snp132,snp138,ljb26_all,gwasCatalog,targetScanS,cosmic70,caddgt20,nci60,clinvar_20140929 -operation g,g,r,r,f,f,f,f,f,f,r,r,f,f,f,f -nastring . -csvout


mv refgene.hg19_multianno.csv $(echo refgene.hg19_multianno.csv | sed "s/\./".$nom_fichier"\./")
mv output.vcf $(echo output.vcf | sed "s/\./".$nom_fichier"\./")
mv output.vcf.idx $(echo output.vcf.idx | sed "s/\./".$nom_fichier"\./")

rm outav
UPTIME14=`uptime`






##log

echo "---------------------------------"
echo "log of the pipeline preGATK"
echo "analyse du ficher = $nom_fichier"
echo "--------------------"
cal
echo "--------------------"
echo "step1"
echo "Uptime = $UPTIME1"
echo "-------------"
echo "step2"
echo "Uptime = $UPTIME2"
echo "------------"
echo "step3"
echo "Uptime = $UPTIME3"
echo "------------"
echo "step4"
echo "Uptime = $UPTIME4"
echo "------------"
echo "step5"
echo "Uptime = $UPTIME5"
echo "------------"
echo "step6"
echo "Uptime = $UPTIME6"
echo "------------"
echo "step7"
echo "Uptime = $UPTIME7"
echo "------------"
echo "step8"
echo "Uptime = $UPTIME8"
echo "------------"
echo "step9"
echo "Uptime = $UPTIME9"
echo "------------"
echo "step10"
echo "Uptime = $UPTIME10"
echo "------------"
echo "step11"
echo "Uptime = $UPTIME11"
echo "------------"
echo "step12"
echo "Uptime = $UPTIME12"
echo "------------"
echo "---------------------------------"
echo "log of the pipeline in ANNOVAR"
echo "--------------------"
echo "convertion ok"
echo "Uptime = $UPTIME13"
echo "------------"
echo "annotation ok"
echo "Uptime = $UPTIME14"
echo "------------"

exit 0

