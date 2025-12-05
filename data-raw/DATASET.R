library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# RF model of BL
RFmodel_BL <- system.file(
  "extdata",
  "RFmodel_BL.rds",
  package = "GAMBLR.predict") %>%
  readRDS

# RF model of FL
RFmodel_FL <- system.file(
  "extdata",
  "RFmodel_FL.rds",
  package = "GAMBLR.predict") %>%
  readRDS

usethis::use_data(RFmodel_BL, overwrite = TRUE)
usethis::use_data(RFmodel_FL, overwrite = TRUE)

# Features of Chapuy classifier with cluster weights
chapuy_features <- list()

chapuy_features$feature_weights <- system.file(
  "extdata",
  "chapuy_weights.tsv",
  package="GAMBLR.predict") %>%
  read_tsv

chapuy_features$ssm_features <- chapuy_features$feature_weights$Feature[
    !str_detect(
        chapuy_features$feature_weights$Feature,
        "SV|AMP|DEL"
    )
]

chapuy_features$cnv_features <- chapuy_features$feature_weights$Feature[
    str_detect(
        chapuy_features$feature_weights$Feature,
        "AMP|DEL"
    )
]

chapuy_features$cnv_features_cytoband <-
chapuy_features$cnv_features[!grepl("p:|q:", chapuy_features$cnv_features)] %>%
as.data.frame %>%
separate(
  .,
  `.`,
  sep = ":",
  into = c("cytoband", "CNV")
)

chapuy_features$cnv_features_arm <-
chapuy_features$cnv_features[grepl("p:|q:", chapuy_features$cnv_features)] %>%
as.data.frame %>%
separate(
  .,
  `.`,
  sep = ":",
  into = c("arm", "CNV")
)

chapuy_features$sv_features <-
chapuy_features$feature_weights$Feature[
    str_detect(
        chapuy_features$feature_weights$Feature,
        "SV"
    )
] %>%
gsub(
    "SV:",
    "",
    .
) %>%
gsub(
    '[/][A-Z0-9]*',
    "",
    .
)

usethis::use_data(chapuy_features, overwrite = TRUE)


# Lacy classifier

# RF model
RFmodel_Lacy <- system.file(
  "extdata",
  "Lacy_rf_model.rds",
  package="GAMBLR.predict") %>%
  readRDS

# Features
lacy_features <- list()

lacy_features$all <- system.file(
  "extdata",
  "lacy_weights.tsv",
  package="GAMBLR.predict") %>%
  read_tsv

lacy_features$cnv <-
lacy_features$all$Feature[
    str_detect(lacy_features$all$Feature, "amp|del")
  ] %>%
  as.data.frame %>%
  dplyr::rename(
    "Gene" = "."
  ) %>%
  dplyr::mutate(
    Feature = Gene,
    CNV = ifelse(
      grepl(
        "amp",
        Gene
      ),
      "AMP",
      "DEL"
    ),
    Dual = ifelse(
      grepl(
        "_OR_",
        Gene
      ),
      "TRUE",
      "FALSE"
    ),
    Gene = gsub(
    '[_][A-Za-z]*',
    "",
    Gene)
  )

lacy_features$shm <-
  lacy_features$all$Feature[
    str_detect(lacy_features$all$Feature, "_S")
  ] %>%
  as.data.frame %>%
  mutate(
    Gene = gsub(
      "_S",
      "",
      `.`
    )
  ) %>%
  `names<-`(c(
    "Feature",
    "Gene"
  )) %>%
  mutate(
    genome_build = "grch37"
  )

lacy_features$grch37_shm <-
lacy_features$shm %>%
    left_join(
        .,
        GAMBLR::grch37_gene_coordinates,
        by=c("Gene"="gene_name")
    ) %>%
    dplyr::select(chromosome, start, end, everything()) %>%
    arrange(chromosome, start, end) %>%
    dplyr::mutate(
      feature_start = start - 2000,
      feature_end = end,
    )

lacy_features$hg38_shm <-
GAMBLR::hg38_gene_coordinates %>%
    dplyr::filter(
        ensembl_gene_id %in% c(
            lacy_features$grch37_shm$ensembl_gene_id,
            "ENSG00000278677",
            "ENSG00000273802",
            "ENSG00000286522",
            "ENSG00000273983")
    ) %>%
    dplyr::mutate(
        Feature = lacy_features$grch37_shm$Feature,
        genome_build = "hg38",
        Gene = gene_name,
        feature_start = start-2000,
        feature_end = end,
    ) %>%
    select(
        colnames(lacy_features$grch37_shm)
    )

lacy_features$hotspots <-
  lacy_features$all$Feature[
    str_detect(lacy_features$all$Feature, "_noncan")
  ] %>%
  gsub(
    "_noncan",
    "",
    .
  )

lacy_features$ssm <-
  lacy_features$all$Feature[
    str_detect(
      lacy_features$all$Feature,
      "_noncan|_S|amp|del",
      negate = TRUE)
  ] %>%
  gsub(
    '[_][0-9]*',
    "",
    .
  )

lacy_features$ssm <-
c(
  lacy_features$ssm,
  lacy_features$cnv %>%
      dplyr::filter(
        Dual == "TRUE"
      ) %>%
      pull(
        Gene
      )
) %>% sort

usethis::use_data(RFmodel_Lacy, overwrite = TRUE)
usethis::use_data(lacy_features, overwrite = TRUE)


# Introduce LymphGenerator features
lymphgenerator_features <- list()

lymphgenerator_features$CNV <- system.file(
  "extdata",
  "lymphgenerator_cnv.tsv",
  package="GAMBLR.predict") %>%
  read_tsv %>%
  arrange(genome_build, chromosome, start)

lymphgenerator_features$SV <- system.file(
  "extdata",
  "lymphgenerator_sv.tsv",
  package="GAMBLR.predict") %>%
  read_tsv %>%
  pull(Tx)

lymphgenerator_features$SV <- gsub(
    "_SV",
    "",
    lymphgenerator_features$SV
)

lymphgenerator_features$SV <- c(
    lymphgenerator_features$SV,
    "MYC"
)

lymphgenerator_features$hotspots <- system.file(
  "extdata",
  "lymphgenerator_hotspots.tsv",
  package="GAMBLR.predict") %>%
  read_tsv %>%
  pull(HOTSPOT)

lymphgenerator_features$SSM <- system.file(
  "extdata",
  "lymphgenerator_ssm.tsv",
  package="GAMBLR.predict") %>%
  read_tsv %>%
  pull(SSM) %>%
  sort

lymphgenerator_features$aSHM <- system.file(
  "extdata",
  "lymphgenerator_ashm.tsv",
  package="GAMBLR.predict") %>%
  read_tsv %>%
  pull(aSHM) %>%
  sort

usethis::use_data(lymphgenerator_features, overwrite = TRUE)


dlbclass_features_all  = construct_reduced_winning_version(
  mutations = system.file("extdata","DLBCL.699.fullGSM.Sep_23_2022.tsv",package="GAMBLR.predict"),
  fisher_test_result = system.file("extdata","fisher_exact_5x2.Sep_23_2022.combined.tsv",package="GAMBLR.predict"))
dlbclass_features = list()
dlbclass_features$SSM = dlbclass_features_all$feature_names$ssm_features

dlbclass_features$SV = dlbclass_features_all$feature_names$sv_features

dlbclass_features$CNV = dlbclass_features_all$feature_names$cnv_features

dlbclass_features$hotspots = dlbclass_features_all$feature_names$hotspot_features

usethis::use_data(dlbclass_features, overwrite = TRUE)

dlbclass_features_all <- read_tsv(
    system.file(
        "extdata",
        "Chapuy_feature_cluster_assoc.tsv",
        package="GAMBLR.predict"
    )
)
usethis::use_data(dlbclass_features_all, overwrite = TRUE)
