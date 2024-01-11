hg19_build_flavours <- c("hg19", "grch37", "hs37d5", "GRCh37")
hg38_build_flavours <- c("hg38", "grch38", "GRCh38")

# Needed for the feature flattening
bcl2_feat <- c("BCL2-intronic", "BCL2-TSS", "BCL2")
myc_feat <- c("MYC-TSS", "MYC", "MYC_IGH")
bcl6_feat <- c("BCL6_SV", "BCL6")
pim1_feat <- c("PIM1", "PIM1-TSS")
btg2_feat <- c("BTG2", "BTG2-intronic")
dtx1_feat <- c("DTX1", "DTX1-TSS")
actb_feat <- c("ACTB", "ACTB-TSS")

features <- list(
	bcl2_feat, myc_feat, bcl6_feat, pim1_feat, btg2_feat, dtx1_feat, actb_feat
)
names <- c(
	"BCL2_any", "MYC_aSHM",
	"BCL6_SV_or_Mut", "PIM1_any",
	"BTG2_any", "DTX1_any",
	"ACTB_any"
)
