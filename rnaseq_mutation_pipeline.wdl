workflow rnaSeqMutationPipeline
	File inputSamplesFile
  	Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
  	File refFasta
  	File refIndex
  	File refDict

  	File gatk
  	File STAR # path to STAR executable
  	String picard

	scatter (sample in inputSamples) {
		sub(input_file, "\\.bam$", ".index")
		if(xxx)
			gunzip -c
		call starPass1 {
			input: GATK=gatk, 
			readFilesCommand=readFilesCommand,
			RefFasta=refFasta, 
			RefIndex=refIndex,  
			sampleName=sample[0],
			left=sample[1], 
			right=sample[2]
		}
	}
	call GenotypeGVCFs {
		input: GATK=gatk, 
		RefFasta=refFasta, 
		RefIndex=refIndex, 
		RefDict=refDict, 
		sampleName="CEUtrio", 
		GVCFs=HaplotypeCallerERC.GVCF
	}
}


STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles hg19.fa \
    --sjdbFileChrStartEnd /path/to/1pass/SJ.out.tab --sjdbOverhang 75 --runThreadN <n>

task starPass1 {
	File left
   	File right
 	String sampleName
 	String readFilesCommand
    File genomeDir # path to genome indices directory 
  	
  	String twopassMode
  	Int memory
    Int diskSpace
    Int threads
    Int preempt

	command {
		${STAR} \
		--genomeDir $(genomeDir} \
		--runThreadN $(threads) \
		--readFilesIn $(left) $(right) \
		${"--readFilesCommand " + readFilesCommand} \
		--outSAMtype BAM SortedByCoordinate \
		--twopassMode Basic \
		--limitBAMsortRAM ${memory} \
		--outFileNamePrefix ${sampleName} \
	}

   output {
        File bam_file = "$(sampleName).Aligned.sortedByCoord.out.bam"
        File bam_index = "$(sampleName).Aligned.sortedByCoord.out.bam.bai"
        File transcriptome_bam = "$(sampleName).Aligned.toTranscriptome.out.bam"
        File chimeric_junctions = "$(sampleName).Chimeric.out.junction"
        File chimeric_bam_file = "$(sampleName).Chimeric.out.sorted.bam"
        File chimeric_bam_index = "$(sampleName).Chimeric.out.sorted.bam.bai"
        File read_counts = "$(sampleName).ReadsPerGene.out.tab"
        File junctions = "$(sampleName).SJ.out.tab"
        File junctions_pass1 = "$(sampleName)._STARpass1/SJ.out.tab"
        Array[File] logs = ["$(sampleName).Log.final.out", "$(sampleName).Log.out", "$(sampleName).Log.progress.out"]
    }

    runtime {
        docker: "broadinstitute/gtex_rnaseq:V8"
        memory: "${memory}GB"
        disks: "local-disk ${diskSpace} HDD"
        cpu: "${threads}"
        preemptible: "${preempt}"
    }
}

task AddOrReplaceReadGroups {
	String sampleName
	File inputBam
	String platform
  
	command {
	    ${picard} AddOrReplaceReadGroups \
	    I=$(inputBam) \
	    O=$(sampleName)_read_groups.bam \
	    SO=coordinate \
	    RGID=id  \
	    RGLB=library \ 
	    RGPL=$(platform) \
	    RGPU=machine \
	    RGSM=$(sampleName) \    
	}

	output {
    	File bamFile = "$(sampleName)_read_groups.bam"
	}
}

task MarkDuplicates {
	File inputBam
	File sampleName

	command {
		$(picard) MarkDuplicates \ 
    	I=$(inputBam) \ 
    	O=$(sampleName)_dedup_reads.bam \
    	M=$(sampleName)_metrics.txt \
    	CREATE_INDEX=true
	}

}

task SplitNCigarReads {

	command {
		$(gatk) -T SplitNCigarReads \
		-R ref.fasta \
		-I dedupped.bam \
		-o split.bam \
		-rf ReassignOneMappingQuality \
		-RMQF 255 \
		-RMQT 60 \
		-U ALLOW_N_CIGAR_READS
	}

}

[Thu Oct 19 06:33:25 EDT 2017] picard.sam.MarkDuplicates INPUT=[/broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/sorted.bam] OUTPUT=/broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/dedupped.bam METRICS_FILE=/broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/mark_duplicates_qc_metrics.txt CREATE_INDEX=true    PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates REMOVE_DUPLICATES=false ASSUME_SORTED=false MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 READ_NAME_REGEX=[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).* OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_MD5_FILE=false


 Program Args: -T SplitNCigarReads -R /seq/regev_genome_portal/RESOURCES/human/Hg19/Hg19.fa -I /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/dedupped.bam -o /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 


 Program Args: -T RealignerTargetCreator -R /seq/regev_genome_portal/RESOURCES/human/Hg19/Hg19.fa -I /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/split.bam --out /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/forIndelRealigner.intervals --known /seq/regev_genome_portal/RESOURCES/human/Hg19/dbsnp_138.b37_Hg19.vcf 

 NFO  07:33:19,553 HelpFormatter - Program Args: -T IndelRealigner -R /seq/regev_genome_portal/RESOURCES/human/Hg19/Hg19.fa -I /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/split.bam -targetIntervals /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/forIndelRealigner.intervals --out /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/realigned.bam -known /seq/regev_genome_portal/RESOURCES/human/Hg19/dbsnp_138.b37_Hg19.vcf 


 INFO  07:48:07,188 HelpFormatter - Program Args: -T BaseRecalibrator -I /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/realigned.bam -R /seq/regev_genome_portal/RESOURCES/human/Hg19/Hg19.fa --out /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/recal_table.table -knownSites /seq/regev_genome_portal/RESOURCES/human/Hg19/dbsnp_138.b37_Hg19.vcf 

INFO  08:01:42,492 HelpFormatter - Program Args: -R /seq/regev_genome_portal/RESOURCES/human/Hg19/Hg19.fa -T PrintReads --out /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/recalibrated.bam -I /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/realigned.bam --BQSR /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/recal_table.table 
INFO  08:01:42,498 HelpFormatter - Executing as kcopipe@uger-c073.broadinstitute.org on Linux 2.6.32-696.6.3.el6.x86_64 amd64; Java 

INFO  08:16:50,334 HelpFormatter - Program Args: -T HaplotypeCaller -R /seq/regev_genome_portal/RESOURCES/human/Hg19/Hg19.fa -I /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/misc/recalibrated.bam -recoverDanglingHeads -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 --out /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/variants.vcf 
NFO  09:16:38,316 HelpFormatter - Program Args: -T VariantFiltration -R /seq/regev_genome_portal/RESOURCES/human/Hg19/Hg19.fa -V /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/variants.vcf -window 35 -cluster 3 -filterName FS -filter FS > 30.0 -filterName QD -filter QD < 2.0 --out /broad/hptmp/jgould/samples/MCF7/RNASEQ_MUTATION_PREMADE/variants_initial_filtering.vcf 
