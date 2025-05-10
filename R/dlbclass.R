
#' Construct reduced 21-dimension feature vector for DLBCLass
#'
#' This function is an R port of some of the pre-processing code
#' in DLBCLass https://github.com/getzlab/DLBCL-Classifier
#' 
#' @param data Data frame of mutation status from DLBCLass supplement
#' @param fisher_test_result Fisher's test result from DLBCLass github repository
#' @param add_missing_features Set to TRUE to fill in missing features with zeroes
#' @returns a list of data frames with the full mutation features, collapsed
#' 21-dimension features and mutation-only (no CNV) features
#' @export
#'
#' @examples
#' 
#' original_dlbclass_features = construct_reduced_winning_version()
#' 
construct_reduced_winning_version <- function(mutations_file = "inst/extdata/DLBCL.699.fullGSM.Sep_23_2022.tsv",
                                              fisher_test_result_file = "inst/extdata/fisher_exact_5x2.Sep_23_2022.combined.tsv",
                                              include_cn = TRUE,
                                              mutation_data,
                                              add_missing_features = FALSE) {



  if(missing(mutations_file)){
    if(missing(mutation_data)){
      stop("Please provide either a mutations file or a mutation_data data frame.")
    }
  }else{
    mutation_data = read.table(mutations,sep="\t",row.names = 1,header=1) 
  }
  
  
  # Transpose data if 'MYD88' is in row names
  if ("MYD88" %in% rownames(mutation_data)) {
    mutation_data <- t(mutation_data) %>% as.data.frame() %>% column_to_rownames("sample_id")
  }
  mutation_data = mutation_data %>% dplyr::select(-any_of(c("PLOIDY","PURITY","COO")),-ends_with("CCF"))
  #print(colnames(mutation_data))
  # Read the q-value data
  qval_df <- read.csv(fisher_test_result_file, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
  
  # Drop genes not in qval_df
  genes_to_drop <- setdiff(colnames(mutation_data), rownames(qval_df))
  #mutation_data = dplyr::select(mutation_data, -any_of(genes_to_drop))
  
  # Add missing features if required
  if (add_missing_features) {
    missing_features <- rownames(qval_df)[!(rownames(qval_df) %in% colnames(mutation_data)) & qval_df$q <= 0.10]
    for (feature in missing_features) {
      cat("Feature", feature, "not found. Adding with full zeros.\n")
      mutation_data[, feature] <- 0
    }
  }
  

  saved_rownames = rownames(mutation_data)
  mutation_data <- as_tibble(mutation_data)
  
  # Convert all columns to integers
  mutation_data <- mutation_data %>%
    mutate(across(everything(), ~ as.integer(as.numeric(.)))) %>%
    as.data.frame()
  rownames(mutation_data) = saved_rownames
  full_data = mutation_data
  full_data[is.na(full_data)] = 0
  data_no_cn = full_data %>% dplyr::select(-ends_with(".AMP"),-ends_with(".DEL"))
  #drop non-mutation columns
  
  full_data = full_data %>% dplyr::select(-ends_with("CCF"))
  #mutation_data <- mutation_data[, !colnames(mutation_data) %in% genes_to_drop]

  # Aggregate specific features into vectors
  print(grep("BCL6",colnames(mutation_data),value=T))
  BCL6_ALT <- rowSums(select(mutation_data, any_of(c("BCL6_SV","SV.BCL6", "BCL6"))), na.rm = TRUE)
  NOTCH2_vec <- rowSums(select(mutation_data, any_of(c("NOTCH2","NOTCH2HOTSPOT", "SPEN", "DTX1"))), na.rm = TRUE)
  M88O_vec <- rowSums(select(mutation_data, any_of(c("MYD88","MYD88.OTHER", "TNFAIP3", "TNIP1", "BCL10", "NFKBIE"))), na.rm = TRUE)
  
  CD70_vec <- rowSums(select(mutation_data, any_of(c("CD70", "FAS", "CD58", "B2M", "FADD", "HLA.B","HLA-B"))), na.rm = TRUE)
  

  
  
  # Additional vectors for other clusters
  #BCL2_combined <- rowSums(mutation_data[, c("BCL2", "SV.BCL2")], na.rm = TRUE)
  BCL2_combined <- rowSums(select(mutation_data, (starts_with("BCL2"))), na.rm = TRUE)

  #CREBBP_vec <- rowSums(mutation_data[, c("CREBBP", "EZH2", "KMT2D", "EP300")], na.rm = TRUE)
  CREBBP_vec <- rowSums(select(mutation_data, starts_with("CREBBP"),
                                                   starts_with("EZH2"),
                                                   starts_with("KMT2D"),
                                                   starts_with("EP300")), na.rm = TRUE)

  if(include_cn){
    C1_vec4 <- rowSums(select(mutation_data, any_of(c("UBE2A", "TMEM30A", "ZEB2", "GNAI2", "X5P.AMP", "POU2F2", "IKZF3", "X3Q28.DEL", 
                              "EBF1", "LYN", "BCL7A", "CXCR4", "CCDC27", "TUBGCP5", "SMG7", "RHOA", "BTG2"))), na.rm = TRUE)
        TP53_biallelic <- rowSums(mutation_data[, c("TP53", "X17P.DEL")], na.rm = TRUE)
       X21Q_AMP <- mutation_data[, "X21Q.AMP"]
      GNA13_vec <- rowSums(mutation_data[, c("GNA13", "TNFRSF14", "MAP2K1", "MEF2B", "IRF8", "HVCN1", 
                                "GNAI2", "MEF2C", "SOCS1", "EEF1A1", "RAC2", "X12Q.AMP", 
                                "POU2AF1", "X6Q14.1.DEL")], na.rm = TRUE)
      Sum_C2_ARM <- rowSums(mutation_data[, c("X17P.DEL", "X21Q.AMP", "X11Q.AMP", "X6P.AMP", "X11P.AMP", "X6Q.DEL", "X7P.AMP", "X13Q.AMP", 
                                 "X7Q.AMP", "X3Q.AMP", "X5P.AMP", "X18P.AMP", "X3P.AMP", "X19Q.AMP", "X9Q.AMP", "X12P.AMP", "X12Q.AMP")], na.rm = TRUE)
  
      Sum_C2_FOCAL <- rowSums(mutation_data[, c("X1P36.11.DEL", "X1P31.1.DEL", "X1P13.1.DEL", "X2Q22.2.DEL", "X16Q12.1.DEL", "X14Q32.31.DEL", 
                                   "X1P36.32.DEL", "X15Q15.3.DEL", "X4Q21.22.DEL", "X9P21.3.DEL", "X8Q24.22.AMP", "X12P13.2.DEL", 
                                   "X2P16.1.AMP", "X8Q12.1.DEL", "X19P13.2.DEL", "X17Q25.1.DEL", "X1Q42.12.DEL", "X3P21.31.DEL", 
                                   "X18Q23.DEL", "X19P13.3.DEL", "X13Q34.DEL", "X7Q22.1.AMP", "X10Q23.31.DEL", "X9P24.1.AMP", 
                                   "X3Q28.AMP", "X11Q23.3.AMP", "X17Q24.3.AMP", "X3Q28.DEL", "X13Q14.2.DEL", "X18Q21.32.AMP", 
                                   "X19Q13.32.DEL", "X6P21.1.AMP", "X18Q22.2.AMP", "EP300", "ZNF423", "CD274")], na.rm = TRUE)
      PTEN <- rowSums(mutation_data[, c("PTEN", "X10Q23.31.DEL", "X13Q14.2.DEL")], na.rm = TRUE)
      Sum_C5_CNA <- rowSums(mutation_data[, c("X18Q.AMP", "X3Q.AMP", "X3P.AMP", "X19Q13.42.AMP", "X6Q21.DEL", 
                                 "X18P.AMP", "X19Q.AMP", "X8Q12.1.DEL", "X6Q14.1.DEL", "X19P13.2.DEL", 
                                 "X9P21.3.DEL", "X18Q21.32.AMP", "X18Q22.2.AMP", "X1Q42.12.DEL", "X1Q32.1.AMP", "X6P21.33.DEL")], na.rm = TRUE)
  
    CN_2P16_1_AMP <- mutation_data[, "X2P16.1.AMP"]
  
  }else{
      C1_vec4 <- rowSums(select(mutation_data,
                          starts_with("UBE2A"), 
                          starts_with("TMEM30A"),
                          starts_with("ZEB2"), starts_with("GNAI2"),
                          starts_with("POU2F2"), starts_with("IKZF3"), 
                          starts_with("EBF1"), starts_with("LYN"), 
                          starts_with("BCL7A"), 
                          starts_with("CXCR4"), starts_with("CCDC27"), 
                          starts_with("TUBGCP5"), 
                              starts_with("SMG7"), 
                              starts_with("RHOA"), 
                              starts_with("BTG2")), na.rm = TRUE)
      TP53_biallelic <- rowSums(mutation_data[, c("TP53"),drop=FALSE], na.rm = TRUE)

      PTEN <- rowSums(mutation_data[, c("PTEN"),drop=FALSE], na.rm = TRUE)
      #GNA13_vec <- rowSums(mutation_data[, c("GNA13", "TNFRSF14", "MAP2K1", "MEF2B", "IRF8", "HVCN1", 
      #                          "GNAI2", "MEF2C", "SOCS1", "EEF1A1", "RAC2", 
      #                          "POU2AF1")], na.rm = TRUE)
      GNA13_vec <- rowSums(select(mutation_data, any_of(c("GNA13", "TNFRSF14","TNFRSF14HOTSPOT",
                                "MAP2K1", "MEF2B","MEF2BHOTSPOT", "IRF8","IRF8HOTSPOT", "HVCN1", "HVCN1HOTSPOT",
                                "GNAI2","GNAI2HOTSPOT",
                                "MEF2C","MEF2CHOTSPOT", 
                                "SOCS1","SOCS1HOTSPOT",
                                "EEF1A1","EEF1A1HOTSPOT",
                                "RAC2", "RAC2HOTSPOT",
                                "POU2AF1","POU2AF1HOTSPOT"))), na.rm = TRUE)
  }

  
  
  SV_MYC <- select(mutation_data, any_of(c("SV.MYC", "MYC","MYC_SV"))) %>% rowSums(na.rm = TRUE)
  
  Hist_comp <- rowSums(select(mutation_data, any_of(c("HIST1H2AC","HIST1H2ACHOTSPOT", 
                                                    "HIST1H1E", "HIST1H1EHOTSPOT","HIST1H1B",
                                                     "HIST1H2AM","HIST1H2AMHOTSPOT", "HIST1H1C",
                                                     "HIST1H1CHOTSPOT","HIST1H1D", "HIST1H1DHOTSPOT",
                                                     "HIST1H2BC","HIST1H2BCHOTSPOT"))), na.rm = TRUE)
  SGK1_vec <- rowSums(select(mutation_data, any_of(c("SGK1", "TET2", "SGK1HOTSPOT","NFKBIA","NFKBIAHOTSPOT",
                              "STAT3","STAT3HOTSPOT", "PTPN6", "BRAF","BRAFHOTSPOT", "KRAS", "KRASHOTSPOT",
                               "CD83","CD83HOTSPOT", "SF3B1","SF3B1HOTSPOT", "CD274", "MEF2C","MEF2CHOTSPOT", "KLHL6","KLHL6HOTSPOT", "CXCR4","CXCR4HOTSPOT", "PTEN", 
                               "RAC2", "SESN3", "SOCS1","SOCS1HOTSPOT", "METAP1D"))), na.rm = TRUE)
  #SGK1_vec <- rowSums(mutation_data[, c("SGK1", "TET2", "NFKBIA", "STAT3", "PTPN6", "BRAF", "KRAS", 
  #                             "CD83", "SF3B1", "CD274", "MEF2C", "KLHL6", "CXCR4", "PTEN", 
  #                             "RAC2", "SESN3", "SOCS1", "METAP1D")], na.rm = TRUE)
  #DUSP2_vec <- rowSums(mutation_data[, c("DUSP2", "ZFP36L1", "CRIP1", "ACTB", "LTB", "YY1", "PABPC1")], na.rm = TRUE)
  DUSP2_vec <- rowSums(select(mutation_data, any_of(c("DUSP2", "DUSP2HOTSPOT", "ZFP36L1","ZFP36L1HOTSPOT", "CRIP1", "ACTB", "LTB", "YY1", "PABPC1","PABPC1HOTSPOT"))), na.rm = TRUE)

  #TBL1XR1_vec <- rowSums(mutation_data[, c("TBL1XR1", "PIM1", "PRDM1", "ETV6", "ZC3H12A", "BTG1", "BTG2", 
  #                                "IGLL5", "TMSB4X", "GRHPR", "HLA.C", "MYD88", "TOX", "LYN", 
  #                                "POU2F2", "IKZF3", "HLA.A", "ZFP36L1", "CARD11", "SF3B1", 
  #                                "HLA.B", "IRF2BP2", "OSBPL10", "ATP2A2", "PIM2", "IRF4", "BCL11A", 
  #                                "METAP1D", "ETS1", "CCDC27")], na.rm = TRUE)
    TBL1XR1_vec <- rowSums(select(mutation_data, any_of(c("TBL1XR1", "PIM1", "PIM1HOTSPOT",
                                                          "PRDM1","PRDM1HOTSPOT",
                                                          "ETV6","ETV6HOTSPOT", "ZC3H12A",
                                                          "BTG1", "BTG1HOTSPOT", "BTG2", "BTG2HOTSPOT",
                                  "IGLL5", "IGLL5HOTSPOT",
                                  "TMSB4X","TMSB4XHOTSPOT", 
                                  "GRHPR","GRHPRHOTSPOT","HLA.C","HLA-C","HLA-CHOTSPOT",
                                  "MYD88", "TOX", "LYN", 
                                  "POU2F2","POU2F2HOTSPOT",
                                  "IKZF3","IKZF3HOTSPOT", 
                                  "HLA.A", "ZFP36L1","HLA-A","HLA-AHOTSPOT",
                                  "CARD11", "CARD11HOTSPOT","SF3B1", 
                                  "HLA.B","HLA-B","HLA-BHOTSPOT", "IRF2BP2", "OSBPL10", "ATP2A2", "PIM2", "IRF4", "BCL11A", 
                                  "METAP1D", "ETS1", "CCDC27"))), na.rm = TRUE)

  MYD88_L265P_CD79B <- rowSums(select(mutation_data, any_of(c("MYD88.L265P","MYD88HOTSPOT", "CD79B","CD79BHOTSPOT"))), na.rm = TRUE)
  
  if(include_cn){
    # Combine all vectors into a data frame
  reduced_data <- data.frame(
    BCL6_ALT = BCL6_ALT,
    NOTCH2_vec = NOTCH2_vec,
    M88O_vec = M88O_vec,
    C1_vec4 = C1_vec4,
    CD70_vec = CD70_vec,
    TP53_biallelic = TP53_biallelic,
    X21Q_AMP = X21Q_AMP,
    Sum_C2_ARM = Sum_C2_ARM,
    Sum_C2_FOCAL = Sum_C2_FOCAL,
    BCL2_combined = BCL2_combined,
    CREBBP_vec = CREBBP_vec,
    GNA13_vec = GNA13_vec,
    PTEN = PTEN,
    SV_MYC = SV_MYC,
    Hist_comp = Hist_comp,
    SGK1_vec = SGK1_vec,
    DUSP2_vec = DUSP2_vec,
    CN_2P16_1_AMP = CN_2P16_1_AMP,
    TBL1XR1_vec = TBL1XR1_vec,
    MYD88_L265P_CD79B = MYD88_L265P_CD79B,
    Sum_C5_CNA = Sum_C5_CNA
    )
  }else{
    # Combine all vectors into a data frame
    reduced_data <- data.frame(
      BCL6_ALT = BCL6_ALT,
      NOTCH2_vec = NOTCH2_vec,
      M88O_vec = M88O_vec,
      C1_vec4 = C1_vec4,
      CD70_vec = CD70_vec,
      TP53_biallelic = TP53_biallelic,
      #X21Q_AMP = X21Q_AMP,
      #Sum_C2_ARM = Sum_C2_ARM,
      #Sum_C2_FOCAL = Sum_C2_FOCAL,
      BCL2_combined = BCL2_combined,
      CREBBP_vec = CREBBP_vec,
      GNA13_vec = GNA13_vec,
      PTEN = PTEN,
      SV_MYC = SV_MYC,
      Hist_comp = Hist_comp,
      SGK1_vec = SGK1_vec,
      DUSP2_vec = DUSP2_vec,
      #CN_2P16_1_AMP = CN_2P16_1_AMP,
      TBL1XR1_vec = TBL1XR1_vec,
      MYD88_L265P_CD79B = MYD88_L265P_CD79B
      #Sum_C5_CNA = Sum_C5_CNA
    )
  }
  
  
  ssm = gsub("\\.","-",colnames(full_data))
  ssm = ssm[!grepl("-DEL|-AMP|SV-|MYD88",ssm)]
  hotspot = "MYD88"
  sv = c("BCL2","BCL6","MYC")
  cnv = grep("\\.AMP|\\.DEL",colnames(full_data),value=T)
  cnv = gsub("^X","",cnv)
  cnv = gsub("\\.AMP","-amp",cnv)
  cnv = gsub("\\.DEL","-del",cnv)
  cnv = gsub("P","p",cnv)
  cnv = gsub("Q","q",cnv)
  return(list(reduced=reduced_data,
              full=full_data,
              no_cn=data_no_cn,
              feature_names=list(ssm_features=ssm,hotspot_features=hotspot,cnv_features=cnv,sv_features=sv)))
}
