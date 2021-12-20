#!/usr/bin/env nextflow



// defind parameter


// defind channels


Channel.from([["/home/labs/aharoni/yedidyab/counts_pipeline_training/heli_final/FW-1_S53_L001_R1_001.fastq.gz", "/home/labs/aharoni/yedidyab/counts_pipeline_training/heli_final/FW-1_S53_L002_R1_001.fastq.gz"], ["/home/labs/aharoni/yedidyab/counts_pipeline_training/heli_final/FW-1_S53_L001_R2_001.fastq.gz", "/home/labs/aharoni/yedidyab/counts_pipeline_training/heli_final/FW-1_S53_L002_R2_001.fastq.gz"]],[["/home/labs/aharoni/yedidyab/counts_pipeline_training/heli_final/FW-2_S54_L001_R1_001.fastq.gz", "/home/labs/aharoni/yedidyab/counts_pipeline_training/heli_final/FW-2_S54_L002_R1_001.fastq.gz"], ["/home/labs/aharoni/yedidyab/counts_pipeline_training/heli_final/FW-2_S54_L001_R2_001.fastq.gz", "/home/labs/aharoni/yedidyab/counts_pipeline_training/heli_final/FW-2_S54_L002_R2_001.fastq.gz"]]).into{ fastq_ch; fastq_ch1 }



Channel
    .fromPath(params.gff)
    .into { gff1; gff2 }



println """\


         R N A  -  S E Q  -  P I P E L I N E
         ====================================
         Work directory: ${params.mainDir}
         fasta file    : ${params.fasta}
         gff file      : ${params.gff}
         seq zip file  : ${params.seqZipFile}
         """
         .stripIndent()







/*
 * STEP 2 - merge L1 & L2
 */

process mergeLfiles {
  publishDir "${params.mainDir}/merge", mode: 'copy'

  input:
  tuple path(R1), val(R2) from fastq_ch
  
  output:
  tuple file('*R1.fq.gz'), file('*R2.fq.gz') into fastq
  
  script:
  """
  name=`echo ${R1[0]} |awk -F_ '{print \$1}'`
  cat ${R1[0]} ${R1[1]} > \${name}__R1.fq.gz
  cat ${R2[0]} ${R2[1]} > \${name}__R2.fq.gz
  """


}

/*
 * STEP 3 - transfer the UMIs from R2 to R1
 */

process umitools {

  input:
  tuple path(R1), path(R2) from fastq

  output:
  file "*" into (umi)

  script:
  """
  name=`echo $R1 | awk -F__ '{print \$1}'`
  umi_tools extract --bc-pattern=NNNN --stdin=$R2 --read2-in=$R1 --stdout=\${name}_UMI.fq.gz --read2-stdout
  """
}


/*
 * STEP 4 - trim with trim galore
 */

process trimGalore {

  input:
  file fq from umi

  output:
  file '*trimmed.fq.gz' into (clean_fastq_ch)

  script:
  """
  trim_galore ${fq} --no_report_file
  """
}

/*
 * STEP 4.2 - PolyA with trim galore
 */

process trimGalorePolyA {

  input:
  file x from clean_fastq_ch

  output:
  file '*' into (clean_fastq_ch2)

  script:
  """
  trim_galore --polyA ${x} --no_report_file
  """
}

/*
 * STEP 5 - move Umi to the end
 */

process moveUmi {

  input:
  file x from clean_fastq_ch2

  output:
  file '*_trimmed' into (to_mapping_ch)

  shell:
  '''
  bash !{params.mainDir}/templates/move.sh !{x} > !{x}_trimmed
  '''
  
}

/*
 * STEP 6.1 - STAR DB 
 */
 
process starDb {
  
  output:
  file 'db_file' into starDB_ch
  
  script:
  """
  mkdir db_file
  STAR --runThreadN 20 --genomeDir db_file --runMode genomeGenerate --genomeSAindexNbases 10 --genomeFastaFiles ${params.fasta} 
  """
}

/*
 * STEP 6.2 - STAR mapping
 */
 
process mappingStar {
  publishDir "${params.mainDir}/mapping", mode: 'copy'


  input:
  file fastq from to_mapping_ch
  file starDB from starDB_ch

  output:
  file "*.sortedByCoord.out.bam" into (bam_ch_to_count,bam_ch_to_umi2,bam_test)
  
  script:
  """
  STAR --genomeDir $starDB \
  --runThreadN 20 \
  --readFilesIn $fastq \
  --outFileNamePrefix "$fastq"__ \
  --outSAMtype BAM SortedByCoordinate\
  --outSAMunmapped Within\
  --outSAMattributes Standard\
  --readFilesCommand zcat 
  """
  
}


/*
 * STEP 7 - subread - feature count
 */
 
process subReadBefore {

  input:
  file(gff_file) from gff1
  file(bam_files) from bam_ch_to_count.collect()
  
  output:
  file("*") into (feature_counts_ch)
  
  script:
  """
  featureCounts \
  -a $gff_file \
  -t gene \
  -g ID \
  -o featureCounts.txt ${bam_files}
  """
}

/*
 * STEP 8 - Umitools 2
 */
 
process umitools2 {
  conda '/home/labs/aharoni/yedidyab/anaconda3/envs/umitools'

  input:
  file bam_file from bam_ch_to_umi2

  output:
  file "*.deduped.bam" into (umi2,umi_test)

  script:
  """
  name=`echo $bam_file |awk -F__ '{print \$1}'`
  
  samtools index -@ 2 ${bam_file}
  umi_tools dedup --stdin=${bam_file} --log=\${name}.deduped.log > ${bam_file}.deduped.bam
  
  """
}

/*
 * STEP 9 - count feture after remove duplicte
 */
 
process subReadAfter {

  input:
  file(gff_file) from gff2
  file(bam_files) from umi2.collect()
  
  output:
  file("*") into (feature_counts_ch2)
  
  script:
  """
  featureCounts \
  -a $gff_file \
  -t gene \
  -g ID \
  -o featureCounts.txt ${bam_files}
  """
}

bam_test.subscribe { println "bam_test: $it" }
umi_test.subscribe { println "bam_test: $it" }
