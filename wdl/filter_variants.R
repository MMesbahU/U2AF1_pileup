
### To Run
## Rscript filter_variants.R sample_id.U2AF1pileup.txt.gz u2af1_hotspot_mut.annot_hg38.tsv sample_id.U2AF1_CH_mutation.maxDP.csv sample_id.U2AF1_CH_mutation.avgDP.csv
###

########################### U2AF1 ############
#### Filter CH variants: DP>=20 & AD.Alt>=5 & ADF_Alt>=1 & ADR_Alt >=1 
### and 
#### For same mutation, take the variant with high DP 
#### or take the avg.
## Get User input
ARGs <- commandArgs(TRUE)

pileup_file <- ARGs[1]

annotation_file <- ARGs[2]

max_dp.output_csv <- ARGs[3]

avg_dp.output_csv <- ARGs[4]
  
##
#### Load packages
library(data.table)
library(stringr)
library(dplyr)
#####

## hotspot mutations
annot_u2af1 <- fread(annotation_file, header = T)
##
####### Read Pileup file
pileups <- fread(pileup_file, header = T)

## Merge with Annotation file
pileups_annot <- merge(pileups[,c(1:5,7:10)],
                       annot_u2af1[,c(1,5,6,4,3)], 
                       by="varID")

# Supporting Forward Reads for Ref and Alt allele
pileups_annot$ADF_Ref <- as.numeric(str_split_fixed(string = pileups_annot$ADF, 
                                                    pattern = "[,]",n = 2)[,1])
pileups_annot$ADF_Alt <- as.numeric(str_split_fixed(string = pileups_annot$ADF, 
                                                    pattern = "[,]",n = 2)[,2])

# Supporting Forward Reads for Ref and Alt allele
pileups_annot$ADR_Ref <- as.numeric(str_split_fixed(string = pileups_annot$ADR, 
                                                    pattern = "[,]",n = 2)[,1])
pileups_annot$ADR_Alt <- as.numeric(str_split_fixed(string = pileups_annot$ADR, 
                                                    pattern = "[,]",n = 2)[,2])
  ## Total Reads Supporting ALT allele 
pileups_annot$AD.Alt <- pileups_annot$ADF_Alt + pileups_annot$ADR_Alt

  ## Filter: DP>=20 & AD.Alt>=5 & ADF_Alt>=1 & ADR_Alt >=1
pileups_annot <- subset(pileups_annot, 
                        pileups_annot$DP>=20 & 
                          pileups_annot$AD.Alt>=5 & 
                          pileups_annot$ADF_Alt>=1 & 
                          pileups_annot$ADR_Alt>=1)

  ## Variant Allele Vractions (VAF)
pileups_annot$VAF <- pileups_annot$AD.Alt/pileups_annot$DP

  ## Find duplicates
pileups_annot$Sample_NonsynOI <- paste(pileups_annot$Sample,
                                       pileups_annot$NonsynOI,
                                   sep = "_")

  ## Remove duplicates
pileups_annot <- pileups_annot[order(pileups_annot$NonsynOI,
                                     pileups_annot$DP,
                                     pileups_annot$AD.Alt,
                                      decreasing = T  ),]
  # Keep Higher DP/AD variant
pileups_annot.v1 <- pileups_annot[!duplicated(pileups_annot$Sample_NonsynOI),]

### Save outputs in CSV format
fwrite(pileups_annot.v1, 
	max_dp.output_csv, 
       row.names = F, 
       col.names = T, 
       sep=",", 
       quote = T)
###########################
  ### Alternative options
  ## take average where supporting reads present in both locations
  # take sum and take avg.
samples_u2af1 <- as.data.frame(table(pileups_annot$Sample))

pileups_annot.v2 <- pileups_annot[pileups_annot$Sample %in% samples_u2af1$Var1[samples_u2af1$Freq==2]  , ]
  # "summarise" function for multiple rows 
pileups_annot.v3 <- pileups_annot.v2 %>%
  group_by(Sample_NonsynOI) %>%
  reframe(Sample,varID,CHROM,POS,REF, ALT,
            Gene,NonsynOI,cosmic96_coding,
            AAChange.refGene,
            DP=sum(DP),
            AD.Alt=sum(AD.Alt),
            VAF=sum(AD.Alt)/sum(DP),
            ADF_Ref = sum(ADF_Ref),
            ADR_Ref=sum(ADR_Ref),
            ADF_Alt = sum(ADF_Alt),
            ADR_Alt=sum(ADR_Alt))
  ## Keep the variant annotations for 
  # "U2AF1" region in chr21:43092356-43108170 region
pileups_annot.v3 <- pileups_annot.v3[order(pileups_annot.v3$Sample_NonsynOI,
                                           pileups_annot.v3$POS,
                                           decreasing = T  ),]

pileups_annot.v3 <- pileups_annot.v3[!duplicated(pileups_annot.v3$Sample_NonsynOI),]

### Save outputs in CSV format
fwrite(pileups_annot.v3, 
       avg_dp.output_csv, 
       row.names = F, 
       col.names = T, 
       sep=",", 
       quote = T)
################## 
