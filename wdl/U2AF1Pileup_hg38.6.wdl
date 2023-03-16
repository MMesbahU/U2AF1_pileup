version 1.0

###################################################################################
## Version 2023-01-10
## Contact Md Mesbah Uddin <mdmesbah@gmail.com>
## Natarajan Lab Broad Institute/MGH

###################################################################################

####################################################################################
# Requires 
	# 1. U2AF1 intervals: "chr21:6484023-6499848,chr21:43092356-43108170" 
    # 2. bcftools docker e.g. "gcr.io/ukbb-analyses/bcftools:1.10.2"
    # 3. Reference fasta, fai and dict
    # 4. Bam file
##### 

#bed coverage
workflow GetU2AF1Pileup {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File bam
    File bai
    Boolean is_bam = true
    # String U2AF1_regions = "chr21:6484023-6499848,chr21:43092356-43108170"
    String U2AF1_regions_pls10bp = "chr21:6484013-6499858,chr21:43092346-43108180"
    String bcftools_docker
    Int boot_disk_size = 3
    Int preemptible_tries = 1
  }

        
 call U2AF1Pileup{
      input:
      input_bam = bam,
      input_bai = bai,
      U2AF1_regions_pls10bp = U2AF1_regions_pls10bp,
      ref_fasta = ref_fasta,
      ref_fai = ref_fai,
      ref_dict = ref_dict,
      is_bam = is_bam,
      docker_image = bcftools_docker,
      preemptible_tries = preemptible_tries
  }

  output {
    File U2AF1mPileup = U2AF1Pileup.outputPileup
  }
}


#
task U2AF1Pileup {
  input {
    # Command parameters
    File input_bam
    File input_bai
    String U2AF1_regions_pls10bp
    File ref_fasta
    File ref_fai
    File ref_dict
    Boolean is_bam
    #String sample_name = basename(input_bam, ".bam")
    #     
    String sample_name = if (is_bam) then basename(input_bam, ".bam") else basename(input_bam, ".cram")
      
    # Runtime parameters
    Int addtional_disk_size = 2
    String machine_mem_size = 3
    String docker_image
    Int preemptible_tries
  }

  #adjust disk size
  Float input_size = size(input_bam, "GB")
  Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB") + size(ref_dict, "GB")
  Float output_size = size(input_bam, "GB") * 0.20
  Int disk_size = ceil(input_size + output_size + addtional_disk_size + ref_size)

  #Calls samtools view to do the conversion
  command {
    set -eo pipefail
    
    echo -e "CHROM\tPOS\tvarID\tREF\tALT\tINFO\tADF\tADR\tDP\tSample" > ${sample_name}.U2AF1pileup.txt
    
    
    ## 
    bcftools mpileup \
    -d 100000 \
    -q 0 \
    -Q 0 \
    -a FORMAT/AD,DP,ADF,ADR,SP,QS,SCR \
    -a INFO/AD,ADF,ADR,SCR \
    -r ~{U2AF1_regions_pls10bp} \
    -Ou \
    -f ~{ref_fasta} ~{input_bam} | \
    bcftools norm -f ~{ref_fasta} -m- |  \
    bcftools query -f '%CHROM\t%POS\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO\t[%ADF]\t[%ADR]\t[%DP]\t[%SAMPLE]\n' >> ${sample_name}.U2AF1pileup.txt
    
    ## 
    gzip ${sample_name}.U2AF1pileup.txt
    
  }

  #Run time attributes:
  #Use a docker with samtools. Set this up as a workspace attribute.
  #cpu of one because no multi-threading is required. This is also default, so don't need to specify.
  #disk_size should equal input size + output size + buffer
  runtime {
    docker: docker_image
    memory: machine_mem_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
    # zones: "us-east1-c us-east1-d us-east1-b us-central1-b us-central1-c us-central1-f us-west1-a us-west1-b us-west1-c us-east4-a us-east4-b us-east4-c"
    # zones: "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-east4-a us-east4-b us-east4-c us-west1-a us-west1-b us-west1-c us-west2-a us-west2-b us-west2-c us-west3-a us-west3-b us-west3-c us-west4-a us-west4-b us-west4-c"
  }
    
  #Outputs a BAM with the same sample name
  output {
    File outputPileup = "${sample_name}.U2AF1pileup.txt.gz"
  }
}