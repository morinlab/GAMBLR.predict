lyseq_genes <- sort(c("BCL2_SV","BCL6_SV","MYD88HOTSPOT",
  "TNFRSF14", "SPEN", "ID3", "ARID1A", "RRAGC",
  "BCL10", "CD58", "NOTCH2", "FCGR2B", "FCRLA",
  "FAS", "CCND1", "BIRC3", "ATM", "KMT2D",
  "STAT6", "BTG1", "FOXO1", "B2M", "MAP2K1",
  "IDH2", "CHD2", "CREBBP", "SOCS1", "IL4R",
  "PLCG2", "TP53", "STAT5B", "STAT3", "CD79B",
  "GNA13", "BCL2", "TCF3", "S1PR2", "JAK3",
  "MEF2B", "DNMT3A", "XPO1", "CXCR4", "SF3B1",
  "PTPN1", "XBP1", "EP300", "MYD88", "SETD2",
  "RHOA", "NFKBIZ", "TBL1XR1", "KLHL6", "TET2",
  "PIM1", "CCND3", "TMEM30A", "PRDM1", "SGK1",
  "TNFAIP3", "CARD11", "POT1", "BRAF", "EZH2",
  "UBR5", "MYC", "NOTCH1", "TRAF2", "P2RY8",
  "BTK", "TP73", "NOL9", "SEMA4A", "PRRC2C",
  "BTG2", "ITPKB", "SEC24C", "EDRF1", "WEE1",
  "MPEG1", "ETS1", "DTX1", "SETD1B", "NFKBIA",
  "ZFP36L1", "CIITA", "MBTPS1", "IRF8", "ACTG1",
  "KLHL14", "MED16", "CD70", "JUNB", "KLF2",
  "CD79A", "DYSF", "DUSP2", "BCL2L1", "PRDM15",
  "RFTN1", "OSBPL10", "EIF4A2", "BCL6", "IRF4",
  "FOXC1", "CD83", "H2BC4", "H1-3", "H1-5",
  "HLA-A", "HLA-B", "PRRC2A", "TBCC", "INTS1",
  "ACTB", "PIK3CG", "PRKDC", "TOX", "CDKN2A",
  "GRHPR", "DDX3X", "PIM2", "UBE2A", "ETV6",
  "MS4A1", "CD19", "HNRNPD", "NFKBIE", "TMSB4X"
))

full_status <- read_tsv(file = paste0(here::here(), "/all_full_status.tsv")) %>% column_to_rownames("sample_id")
tier1_genes = lymphoma_genes %>% filter(DLBCL_Tier==1) %>% pull(Gene)
tier1_genes = tier1_genes[tier1_genes %in% colnames(full_status)]
tier1_genes = unique(c("BCL2_SV","BCL6_SV","MYD88HOTSPOT",tier1_genes))

#dlbcl_meta_clean <- read_tsv(file = paste0(here::here(), "/dlbcl_meta_clean.tsv"))
dlbcl_meta_clean <- read_tsv(file = paste0(here::here(), "/dlbcl_meta_with_dlbclass.tsv")) %>%
  mutate(DLBClass = ifelse(Confidence > 0.7,PredictedCluster,"Other"))

# This code needs to be fixed to work anywhere
lacy_df <- readxl::read_excel("~/git/Lyntegrate/data/bloodbld2019003535-suppl2.xlsx", sheet = 1) %>%
  mutate(Gene = gsub("_.+", "", Gene))
all_lacy_genes <- c(unique(lacy_df$Gene), "MYD88HOTSPOT")
all_lacy_genes <- all_lacy_genes[all_lacy_genes %in% colnames(full_status)]


lacy_df <- lacy_df %>% filter(`Included in statistical analysis` == "Yes")

lacy_genes <- lacy_df %>%
  pull(Gene) %>%
  unique()
lacy_genes <- c(lacy_genes, "MYD88HOTSPOT")
lacy_genes <- lacy_genes[lacy_genes %in% colnames(full_status)]
## make demo data



prototypes <- full_status[c(1:17), ]
rownames(prototypes) <- c("CCS_0004","CCS_0005","CCS_0007","CCS_0009","CCS_0010",
  "BN2_1", "BN2_2", "ST2_1", "ST2_2", "EZB_1", "EZB_2", "MCD_1", "MCD_2", "N1_1", "N1_2", "Other_1", "Other_2")
prototypes[] <- 0
maxval = max(full_status)
prototypes["CCS_0004",c("ETS1","BTG2","KMT2D","BTG1","CD83")] <- maxval
prototypes["CCS_0005",c("DTX1","NOL9","BTG2","PIM1","OSBPL10","TBL1XR1","CD83","ACTG1")] <- maxval
prototypes["CCS_0007",c("BCL6_SV","NOTCH2","DTX1","NOL9","KLF2","BTG2","BCL2_SV","KMT2D","SOCS1","PIM1","HLA-A","CD83","ZFP36L1","DUSP2")] <- maxval
prototypes["CCS_0009",c("DTX1","EDRF1")] <- maxval
prototypes["CCS_0010",c("BCL6_SV","DTX1","PIM1","BTG1","CD83","DUSP2","HIST1H1B")] <- maxval
prototypes["BN2_1", c("BCL6_SV", "TNFAIP3", "KLF2")] <- maxval
prototypes["BN2_2", c("NOTCH2", "BCL6", "SPEN", "UBE2A")] <- maxval
prototypes["ST2_1", c("SGK1", "ZFP36L1", "ACTG1", "STAT3", "CD83")] <- maxval
prototypes["ST2_2", c("TET2", "ITPKB", "NFKBIA", "DUSP2", "JUNB")] <- maxval
prototypes["MCD_1", c("MYD88HOTSPOT", "PIM1", "ETV6")] <- maxval
prototypes["MCD_2", c("CD79B", "PIM1", "BTG1", "GRHPR")] <- maxval
prototypes["EZB_1", c("CREBBP", "BCL2_SV", "KMT2D")] <- maxval
prototypes["EZB_2", c("EZH2", "TNFRSF14", "IRF8")] <- maxval
prototypes["N1_1", c("NOTCH1", "ATM")] <- maxval
prototypes["N1_2", c("NOTCH1", "ID3")] <- maxval
prototypes["Other_1", c("CD83", "TP53", "BTK","NFKBIZ", "STAT6")] <- maxval
prototypes["Other_2", c("RRAGC", "P2RY8", "CDKN2A", "MS4A1", "B2M")] <- maxval




impact_genes <- c(
  "MYD88HOTSPOT", "ABL1", "BCL2", "CEBPA", "ETV6", "HGF", "JUN", "MSH2", "PHF6", "RPTOR",
  "SRSF2", "ZRSR2", "ACTG1", "BCL6", "CHEK1", "EZH2", "HIF1A", "KDM5A", "MSH6", "PIGA", "RRAGC", "STAG1", "H1B",
  "AKT1", "BCOR", "CHEK2", "FAM46C", "HIST1H1B", "KDM5C", "MTOR", "PIK3C2G", "RTEL1", "STAG2", "H1-2", "AKT2",
  "BCORL1", "CIC", "FANCA", "HIST1H1C", "KDM6A", "MUTYH", "PIK3C3", "RUNX1", "STAT3", "H1D", "AKT3", "BCR",
  "CIITA", "FANCC", "HIST1H1D", "KDR", "MYC", "PIK3CA", "RUNX1T1", "STAT5A", "H1E", "ALK", "BIRC3", "CRBN",
  "FANCD2", "HIST1H1E", "KEAP1", "MYCL1", "PIK3CG", "SAMHD1", "STAT5B", "H2AC", "ALOX12B", "BLM", "CREBBP",
  "FAS", "HIST1H2AC", "KIT", "MYCN", "PIK3R1", "SDHA", "STAT6", "H2AG", "AMER1", "BRAF", "CRKL", "FAT1", "HIST1H2AG",
  "KMT2A", "MYD88", "PIK3R2", "SDHB", "STK11", "H2AL", "APC", "BRCA1", "CRLF2", "FBXO11", "HIST1H2AL", "KMT2B", "NBN",
  "PIM1", "SDHC", "SUFU", "H2AM", "AR", "BRCA2", "CSF1R", "FBXW7", "HIST1H2AM", "KMT2C", "NCOR1", "PLCG1", "SDHD",
  "SUZ12", "H2BC", "ARAF", "BRD4", "CSF3R", "FGF19", "HIST1H2BC", "KMT2D", "NCOR2", "PLCG2", "SETBP1", "SYK",
  "H2BC5", "ARHGEF28", "BRIP1", "CTCF", "FGF3", "HIST1H2BD", "KRAS", "NCSTN", "PMS2", "SETD1A", "TBL1XR1",
  "H2BG", "ARID1A", "BTG1", "CTNNB1", "FGF4", "HIST1H2BG", "KSR2", "NF1", "PNRC1", "SETD1B", "TBX3", "H2BJ",
  "ARID1B", "BTK", "CUX1", "FGFR1", "HIST1H2BJ", "LCK", "NF2", "POT1", "SETD2", "TERT", "H2BK", "ARID2", "CALR",
  "CXCR4", "FGFR2", "HIST1H2BK", "LMO1", "NFE2", "PPP2R1A", "SETD3", "TET1", "H2BO", "ARID3A", "CARD11", "CYLD",
  "FGFR3", "HIST1H2BO", "LTB", "NFE2L2", "PRDM1", "SETD4", "TET2", "H3C2", "ARID3B", "CASP8", "DAXX", "FGFR4",
  "HIST1H3B", "MALT1", "NKX2-1", "PRKAR1A", "SETD5", "TET3", "H3C8", "ARID3C", "CBFB", "DDR2", "FLCN", "HIST1H3G",
  "MAP2K1", "NOTCH1", "PTCH1", "SETD6", "TGFBR2", "ARID4A", "CBL", "DDX3X", "FLT1", "HLA-A", "MAP2K2", "NOTCH2",
  "PTEN", "SETD7", "TNFAIP3", "ARID4B", "CCND1", "DIS3", "FLT3", "HNF1A", "MAP2K4", "NOTCH3", "PTPN1", "SETD8",
  "TNFRSF14", "ARID5A", "CCND2", "DNMT3A", "FLT4", "HRAS", "MAP3K1", "NOTCH4", "PTPN11", "SETDB1", "TOP1",
  "ARID5B", "CCND3", "DOT1L", "FOXL2", "ID3", "MAP3K13", "NPM1", "PTPN2", "SETDB2", "TP53", "ASXL1", "CCNE1",
  "DTX1", "FOXO1", "IDH1", "MAP3K14", "NRAS", "RAD21", "SF3B1", "TP63", "ASXL2", "CD274", "DUSP22", "FOXP1",
  "IDH2", "MAPK1", "NSD1", "RAD50", "SGK1", "TRAF2", "ATM", "CD28", "EED", "FURIN", "IGF1", "MAPK3", "NT5C2",
  "RAD51", "SH2B3", "TRAF3", "ATP6AP1", "CD58", "EGFR", "FYN", "IGF1R", "MCL1", "NTRK1", "RAD51B", "SMAD2",
  "TRAF5", "ATP6V1B2", "CD79A", "EGR1", "GATA1", "IGF2", "MDM2", "NTRK2", "RAD51C", "SMAD4", "TSC1", "ATR",
  "CD79B", "EP300", "GATA2", "IKBKE", "MDM4", "NTRK3", "RAD51D", "SMARCA4", "TSC2", "ATRX", "CDC73", "EP400",
  "GATA3", "IKZF1", "MED12", "P2RY8", "RAD52", "SMARCB1", "TSHR", "ATXN2", "CDH1", "EPHA3", "GNA11", "IKZF3",
  "MEF2B", "PAK7", "RAD54L", "SMARCD1", "TYK2", "AURKA", "CDK12", "EPHA5", "GNA12", "IL7R", "MEN1", "PALB2",
  "RAF1", "SMC1A", "U2AF1", "AURKB", "CDK4", "EPHA7", "GNA13", "INPP4B", "MET", "PARP1", "RARA", "SMC3",
  "U2AF2", "AXIN1", "CDK6", "EPHB1", "GNAQ", "IRF1", "MGA", "PAX5", "RB1", "SMG1", "UBR5", "AXL", "CDK8",
  "ERBB2", "GNAS", "IRF4", "MGAM", "PBRM1", "REL", "SMO", "VAV1", "B2M", "CDKN1B", "ERBB3", "GNB1", "IRF8",
  "MITF", "PCBP1", "RET", "SOCS1", "VAV2", "BACH2", "CDKN2A", "ERBB4", "GRIN2A", "IRS2", "MLH1", "PDCD1",
  "RHOA", "SOX2", "VHL", "BAP1", "CDKN2Ap14ARF", "ERG", "GSK3B", "JAK1", "MOB3B", "PDGFRA", "RICTOR", "SP140",
  "WHSC1", "BARD1", "CDKN2Ap16INK 4A", "ESCO2", "HDAC1", "JAK2", "MPEG1", "PDGFRB", "RNF43", "SPEN", "WT1",
  "BCL10", "CDKN2B", "ESR1", "HDAC4", "JAK3", "MPL", "PDPK1", "ROBO1", "SPOP", "XBP1", "BCL11B", "CDKN2C",
  "ETNK1", "HDAC7", "JARID2", "MRE11A", "PDS5B", "ROS1", "SRC", "XPO1"
)

# This is for the model that does not consider A53 so no TP53 and no CNV features
lymphgen_genes <- c(
  "ACTB", "ACTG1", "BCL10", "BCL2",
  "BCL2L1", "BCL6", "BTG1", "BTG2",
  "CD70", "CD79B", "CD83", "CREBBP",
  "DDX3X", "DTX1", "DUSP2", "EDRF1",
  "EIF4A2", "EP300", "ETS1", "ETV6",
  "EZH2", "FOXC1", "GRHPR", "HIST1H1B",
  "HIST1H2BC", "HLA-A", "HLA-B", "ID3",
  "IRF8", "ITPKB", "JUNB", "KLF2",
  "KLHL14", "KMT2D", "MBTPS1", "MED16",
  "MEF2B", "MPEG1", "MYD88HOTSPOT", "NFKBIA",
  "NOL9", "NOTCH1", "NOTCH2", "OSBPL10",
  "PIM1", "PIM2", "PRDM1", "PRRC2A",
  "PPRC2C", "RFTN1", "SEC24C", "SETD1B",
  "SGK1", "SOCS1", "SPEN", "STAT3",
  "TBL1XR1", "TET2", "TNFAIP3", "TNFRSF14",
  "UBE2A", "WEE1", "ZFP36L1", "BCL2_SV", "BCL6_SV"
)
lyseq_genes <- sort(lyseq_genes[lyseq_genes %in% colnames(full_status)])


panels <- list(
  Coyle = tier1_genes,
  Lacy = lacy_genes,
  Lymphgen = lymphgen_genes[lymphgen_genes %in% colnames(full_status)],
  Lymphplex = c(
    "BCL2", "BCL6", "ARID1A", "B2M", "BTG1", "BTG2", "CCND3",
    "CD70", "CD79B", "CIITA", "CREBBP", "DDX3X", "DTX1",
    "DUSP2", "EP300", "EZH2", "FAS", "GNA13", "IRF4",
    "IRF8", "KMT2D", "MPEG1", "MYD88", "NOTCH1", "NOTCH2",
    "PIM1", "PRDM1", "SGK1", "SOCS1",
    "STAT3", "STAT6", "TBL1XR1", "TET2", "TNFAIP3", "TNFRSF14", "ZFP36L1"
  ),
  Lyseq = lyseq_genes,
  `Lyseq Tier 1` = lyseq_genes[lyseq_genes %in% tier1_genes],
  `Lyseq Tier 1 (no TP53)` = lyseq_genes[lyseq_genes != "TP53"],
  `MSK Impact` = impact_genes[impact_genes %in% colnames(full_status)]
)


demo_samples <- rownames(prototypes)