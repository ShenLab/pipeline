## Resource Directories
export EXOMPPLN="/home/local/users/jw/wes_pipeline/Exome-vc-pipeline-V3-master_AshGit" # Directory containing pipeline shell scripts
export EXOMRES="/home/local/users/jw/resources" # Directory containing resources/references for pipeline

#jar files  and directories for software
GATKJAR="/home/local/users/jw/bin/GenomeAnalysisTK.jar" #Current GATK jar file
PICARD="/home/local/users/jw/bin/picard.jar" #directory containing Picard jar files
SNPEFF="/home/local/users/jw/bin/snpEff.jar" # Current snpEff jar file

## References
export BUILD="b38" # shorthand for build
export REF="$EXOMRES/reference_genomes/GRCh38/Homo_sapiens_assembly38.fasta" 
export HAPMAP="$EXOMRES/references/GRCh38/hapmap_3.3.hg38.vcf.gz" # hapmap vcf from GATK
export INDEL="$EXOMRES/references/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" # Gold standard INDEL reference from GATK
export TGVCF="$EXOMRES/references/GRCh38/1000G_omni2.5.hg38.vcf.gz" 
export DBSNP="$EXOMRES/references/GRCh38/dbsnp_146.hg38.vcf.gz" # dbSNP vcf from GATK
export ONEKG="$EXOMRES/references/GRCh38/1000G_phase1.snps.high_confidence.hg38.vcf.gz" # 1000 genome SNPs vcf
export ANNOVAR='/home/local/users/jw/software_packages/annovar'
export ANNHDB='/home/local/users/jw/software_packages/annovar/humandb' #Location of annovar databases
export HUMANREF="$EXOMRES/human_g1k_v37.fasta" # human 1000 genome assembly from GATK


#Capture Kit Target Files
export TGTCODES="AgtV2:AgtV4:AgtV5:AgtV5UTR:NbgV2:NbgV3:IllTS:BigTgt:RefSeq:VCRv2"

