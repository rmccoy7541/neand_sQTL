library(tidyr)
library(dplyr)

introns <- read.table("../intron_counts/GTEx_v8_junctions_nohead.gct.gz",stringsAsFactors=FALSE,header=TRUE)
sqtl <- read.table("../sqtl/sqtl_test.txt",stringsAsFactors=FALSE,header=TRUE)
vcf <- read.table("../vcf/vcf_for_merge.txt",stringsAsFactors=FALSE)

# change intron cluster names to match names in sQTL file
introns2 <- as.data.frame(sapply(introns, gsub, pattern="_", replacement=":"))

# separate sQTL intron cluster ID field to match with IDs from introns file
sqtl2 <- separate(sqtl, phenotype_id, c("cluster_pos","cluster_id"), sep=":clu_", remove=TRUE)

# join introns and sqtl file by intron cluster ID
# if there are multiple variants in the sqtl file that correspond to one intron cluster,
# duplicate the line from the introns file
combined <- inner_join(sqtl2, introns2, by=c("cluster_pos"="Name"))

write.table(combined, file="combined.txt", sep="\t", quote=FALSE)

# join intron-sqtl and vcf dataframes by variant ID
# if there are multiple intron clusters that correspond to one variant,
# duplicate the line from the vcf
combined2 <- inner_join(combined, vcf, by=c("variant_id"="V1"))

write.table(combined2, file="combined2.txt", sep="\t", quote=FALSE)

# combine identical sample columns from vcf and intron files
# want to end up with a single sample column with info: GT;intron read count