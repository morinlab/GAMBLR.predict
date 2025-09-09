

#' Run the DLBCLone Shiny App
#'
#' Launches the Shiny application for interactive exploration.
#'
#' @import shiny
#' @import readr
#' @import dplyr
#' @import DT
#' @import shinybusy
#' @import plotly
#' @import shinyjs
#' @import ggplot2
#' @import ggside
#' @import purrr
DLBCLone_shiny <- function(...){
    ## Preamble
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

default_panel <- "Coyle"

#full_status <- read_tsv(file = paste0(here::here(), "/all_full_status.tsv")) %>% column_to_rownames("sample_id")
full_status = read_tsv(file=system.file(package = "GAMBLR.predict", "extdata", "all_full_status.tsv")) %>% column_to_rownames("sample_id")


tier1_genes = GAMBLR.data::lymphoma_genes %>% filter(DLBCL_Tier==1) %>% pull(Gene)
tier1_genes = tier1_genes[tier1_genes %in% colnames(full_status)]
tier1_genes = unique(c("BCL2_SV","BCL6_SV","MYD88HOTSPOT",tier1_genes))

dlbcl_meta_clean = read_tsv(file=system.file(package = "GAMBLR.predict", "extdata", "dlbcl_meta_with_dlbclass.tsv"))


dlbcl_meta_clean = dlbcl_meta_clean %>% 
    filter(lymphgen %in% c("EZB","MCD","BN2","ST2","N1","Other")) %>%
    mutate(DLBClass = ifelse(Confidence > 0.7,PredictedCluster,"Other"))



lacy_df = readxl::read_excel(system.file(package = "GAMBLR.predict", "extdata", "bloodbld2019003535-suppl2.xlsx")) %>%
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
rownames(prototypes) <- c("EX_0004","EX_0005","EX_0007","EX_0009","EX_0010",
  "BN2_1", "BN2_2", "ST2_1", "ST2_2", "EZB_1", "EZB_2", "MCD_1", "MCD_2", "N1_1", "N1_2", "Other_1", "Other_2")
prototypes[] <- 0
maxval = max(full_status)
prototypes["EX_0004",c("ETS1","BTG2","KMT2D","BTG1","CD83")] <- maxval
prototypes["EX_0005",c("DTX1","NOL9","BTG2","PIM1","OSBPL10","TBL1XR1","CD83","ACTG1")] <- maxval
prototypes["EX_0007",c("BCL6_SV","NOTCH2","DTX1","NOL9","KLF2","BTG2","BCL2_SV","KMT2D","SOCS1","PIM1","HLA-A","CD83","ZFP36L1","DUSP2")] <- maxval
prototypes["EX_0009",c("DTX1","EDRF1")] <- maxval
prototypes["EX_0010",c("BCL6_SV","DTX1","PIM1","BTG1","CD83","DUSP2","HIST1H1B")] <- maxval
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



load_gene_panels = function(full_status){

impact_genes <- c(
  "MYD88HOTSPOT", "MYD88", "ABL1", "BCL2", "CEBPA", "ETV6", "HGF", "JUN", "MSH2", "PHF6", "RPTOR",
  "SRSF2", "ZRSR2", "ACTG1", "BCL6", "CHEK1", "EZH2", "HIF1A", "KDM5A", "MSH6", "PIGA", "RRAGC", "STAG1", "H1B",
  "AKT1", "BCOR", "CHEK2", "FAM46C", "HIST1H1B", "KDM5C", "MTOR", "PIK3C2G", "RTEL1", "STAG2", "H1-2", "AKT2",
  "BCORL1", "CIC", "FANCA", "HIST1H1C", "KDM6A", "MUTYH", "PIK3C3", "RUNX1", "STAT3", "H1D", "AKT3", "BCR",
  "CIITA", "FANCC", "HIST1H1D", "KDR", "MYC", "PIK3CA", "RUNX1T1", "STAT5A", "H1E", "ALK", "BIRC3", "CRBN",
  "FANCD2", "HIST1H1E", "KEAP1", "MYCL1", "PIK3CG", "SAMHD1", "STAT5B", "H2AC", "ALOX12B", "BLM", "CREBBP",
  "FAS", "HIST1H2AC", "KIT", "MYCN", "PIK3R1", "SDHA", "STAT6", "H2AG", "AMER1", "BRAF", "CRKL", "FAT1", "HIST1H2AG",
  "KMT2A", "PIK3R2", "SDHB", "STK11", "H2AL", "APC", "BRCA1", "CRLF2", "FBXO11", "HIST1H2AL", "KMT2B", "NBN",
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

lymphgen_genes = read_tsv(file=system.file(package = "GAMBLR.predict", "extdata", "DLBCLass_and_Lymphgen_genes.tsv"))  %>%
 filter(Lymphgen_no_CNV==1) %>%
    pull(gene)

lymphgen_genes <- sort(intersect(lymphgen_genes,colnames(full_status)))
lyseq_genes <- sort(lyseq_genes[lyseq_genes %in% colnames(full_status)])

dlbclass_genes = read_tsv(file=system.file(package = "GAMBLR.predict", "extdata", "DLBCLass_and_Lymphgen_genes.tsv"))  %>%
 filter(DLBCLass_no_CNV==1) %>%
    pull(gene)
dlbclass_absent = setdiff(dlbclass_genes,colnames(full_status))
print("missing")
print(paste(dlbclass_absent,collapse=","))
dlbclass_genes = intersect(dlbclass_genes,colnames(full_status))
lymphgen_genes <- sort(intersect(lymphgen_genes,colnames(full_status)))



panels <- list(
  Coyle = tier1_genes,
  `DLBClass (no CNV)` = dlbclass_genes,
  Lacy = lacy_genes,
  Lymphgen = lymphgen_genes,
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
return(panels)
}

panels = load_gene_panels(full_status)

demo_samples <- rownames(prototypes)

#default_panel <- "Coyle"
defaul_panel <- "Lyseq Tier 1 (no TP53)"

#default_umap <- make_and_annotate_umap(
#  df = full_status %>% select(all_of(tier1_genes)),
#  metadata = dlbcl_meta_clean
#)
#write_tsv(default_umap$df,file="inst/extdata/default_umap_df.tsv")
#default_umap = list(df=read_tsv("inst/extdata/default_umap_df.tsv"))
default_umap_df = read_tsv(file=system.file(package = "GAMBLR.predict", "extdata", "default_umap_df.tsv")) %>% column_to_rownames("sample_id")
default_umap = list(df=default_umap_df)
k_low <- 10
k_high <- 10
default_mode = "Lenient"
truth_col <- "lymphgen"
other_lbl <- "Other"

# derive class set from metadata for the chosen paradigm
meta_classes <- dlbcl_meta_clean[[truth_col]] %>% unique() %>% setdiff(NA) %>% as.character()
meta_classes = sort(meta_classes[meta_classes!=other_lbl])

meta_classes = c(meta_classes,other_lbl)

default_knn <- DLBCLone_KNN(full_status %>% select(all_of(tier1_genes)),
  dlbcl_meta_clean,
  min_k = k_low,
  max_k = k_high,
  optimize_for_other = if_else(default_mode=="Lenient",FALSE,TRUE)
)

pred_dir <- file.path(tempdir(), "predictions")
dir.create(pred_dir, showWarnings = FALSE, recursive = TRUE)
addResourcePath("predictions", pred_dir)


    ui <- fluidPage(
    useShinyjs(),
    add_busy_spinner(spin = "fading-circle", color = "#003366", position = "full-page"),
    # tags$div(style = "display:none;", downloadButton("download_predictions", "Download CSV")),


    tags$div(
        id = "loading-text", style = "color:#003366; font-weight:bold; display:none;",
        "Running model... please wait."
    ),
    tags$script(HTML("
    $(document).on('click', '.download_btn', function() {
        Shiny.setInputValue('download_run_id', this.id, {priority: 'event'});
    });

    Shiny.addCustomMessageHandler('triggerDownload', function(message) {
        Shiny.setInputValue('download_run_id_internal', message.run_id, {priority: 'event'});
        Shiny.setInputValue('go_now', Math.random());  // force a re-trigger
    });
    ")),

    # Application title
    titlePanel("DLBCLone"),
    sidebarLayout(
        sidebarPanel(
        actionButton("visualize", "Update UMAP", icon("refresh")),
        helpText("Generate a quick visualization of samples using the selected feature space"),
        radioButtons(
            "truth_column", "Classification paradigm",
            choices  = c("LymphGen" = "lymphgen", "DLBClass" = "DLBClass"),
            selected = "lymphgen", inline = TRUE
            ),
        radioButtons("mode","Classification Mode",choices=c("Stringent","Lenient"),selected="Lenient",inline=TRUE),
        selectInput("panel",
            label = "Pre-defined gene panel",
            choices = names(panels),
            selected = default_panel
        ),


        radioButtons(
            "predict_source", "Prediction input",
            choices = c("Demo sample", "Custom sample"),
            selected = "Demo sample", inline = TRUE
        ),
        # Old demo selector lives under Demo
    conditionalPanel(
        condition = "input.predict_source == 'Demo sample'",
        selectInput(
        "predict_sample",
        label = "Predict class for a simulated sample",
        choices = demo_samples,
        selected = demo_samples[1]
        )
    ),

    # NEW: Custom sample inputs
    conditionalPanel(
        condition = "input.predict_source == 'Custom sample'",
        textInput("custom_sample_id", "Sample ID", placeholder = "e.g., MY_CASE_001"),
        helpText("Pick features present in your sample (unselected default to 0)."),
        selectizeInput(
        "custom_selected",
        label = "Features present (0/1 style)",
        choices = sort(colnames(full_status)),
        selected = NULL,
        multiple = TRUE,
        options = list(placeholder = 'Type to search features…')
        ),
        numericInput(
        "custom_present_value",
        "Value to assign to present features",
        value = maxval, min = 0, step = 1, max=maxval
        )
    ),

        actionButton("predict", "Predict!", icon("refresh")),
        helpText("Use the current model to predict the class for this sample."),
        helpText("You will need to regenerate the classifier if any settings have been changed."),
        checkboxGroupInput("features",
            label = "Features", inline = T,
            choices = sort(colnames(full_status)),
            selected = tier1_genes
        ),
        checkboxInput("use_core", "Specify a set of core features (higher importance)", value = FALSE),

    conditionalPanel(
    condition = "input.use_core",
    helpText("Core features will be more heavily weighted in the analysis."),
    selectizeInput(
        "core_features",
        label = "Core features",
        choices = sort(colnames(full_status)),
        selected = NULL,
        multiple = TRUE,
        options = list(placeholder = 'Type to search features…', create=FALSE)
    )
    ),
        sliderInput("k_range", label = "Range of K values to try", min = 5, max = 50, step = 5, value = c(k_low, k_high)),
        actionButton("regenerate", "Update Classifier", icon("refresh")),
        actionButton("share_link", "Get permalink for this configuration", icon("link"))

        ),

        mainPanel(
        tabsetPanel(
            tabPanel("UMAP", plotOutput("umap_scatterplot",height = "800px",width= "700px")),
            tabPanel("UMAP (Lymphgen)", plotlyOutput("DLBCLone_KNN_plot_truth")),
            tabPanel("UMAP (DLBCLone)", plotlyOutput("DLBCLone_KNN_plot_prediction")),
            tabPanel("Results overview", plotOutput("alluvial", height = "800px", width = "700px")),
            tabPanel("All DLBCLone assignments", downloadButton("downloadData", "Download"), DTOutput("predictions")),
            tabPanel("Model Log", DTOutput("run_log_table")),
            tabPanel("Sample Neighbors", plotOutput("heatmap",height="500px", width="800px")),
            tabPanel("Prediction Log", DTOutput("sample_log_table"))
        )
        )
    )
    )

    server <- function(input, output, session) {
    dlbclone_result <- reactiveVal(isolate(default_knn))
    dlbclone_pred_result <-reactiveVal()
    umap_result <- reactiveVal(isolate(default_umap))
    run_log <- reactiveVal()
    prediction_store <- reactiveVal(list())
    sample_log <- reactiveVal()
    # Track the most recently predicted sample ID (demo or custom)
    current_predict_id <- reactiveVal(NULL)
    restoring_from_query <- reactiveVal(FALSE)
    pending_core <- reactiveVal(NULL)


    # Which truth column / classes / outgroup label does the current model use?
    current_truth_column <- reactive({
        res <- dlbclone_result()
        res$truth_column %||% input$truth_column %||% "lymphgen"
    })

    current_classes <- reactive({
        res <- dlbclone_result()
        if (!is.null(res$truth_classes)) as.character(res$truth_classes)
        else {
        # fallback: derive from metadata
        vals <- dlbcl_meta_clean[[ current_truth_column() ]] %>% unique() %>% setdiff(NA)
        as.character(vals)
        }
    })

    current_other_label <- reactive({
        # backend normalizes the outgroup score column to 'Other_score'
        # but predicted label value can be any string; default to "Other"
        "Other"
    })
    # Selected core features from UI (distinct from input id to avoid shadowing)
    selected_core <- reactive({
        if (isTRUE(input$use_core) && length(input$core_features)) {
        input$core_features
        } else {
        character()
        }
    })

    `%||%` <- function(x, y) if (is.null(x)) y else x

    # What settings were last used to train the current model?
    committed <- reactiveValues(
        features = colnames(default_knn$features_df),
        core     = character(),
        use_core = FALSE,
        mode     = "Lenient",
        k_range  = c(k_low, k_high),
        panel    = default_panel
    )

    # Build a query string from a named list (handles URL encoding)
    build_query <- function(params) {
    paste(
        vapply(names(params), function(nm) {
        paste0(nm, "=", URLencode(as.character(params[[nm]]), reserved = TRUE))
        }, character(1)),
        collapse = "&"
    )
    }
    # Testing permalink
    # Build absolute base URL to this app
    base_url <- reactive({
        proto <- session$clientData$url_protocol   # e.g. "http:"
        host  <- session$clientData$url_hostname   # e.g. "127.0.0.1"
        port  <- session$clientData$url_port       # e.g. "1234" (may be "")
        path  <- session$clientData$url_pathname   # e.g. "/"
        paste0(proto, "//", host, if (nzchar(port)) paste0(":", port) else "", path)
    })


    init_done <- reactiveVal(FALSE)

    # Is the UI "dirty" vs the committed settings?
    dirty <- reactive({
        if (!isTRUE(init_done())) return(FALSE)  # not dirty until we snapshot once

        cur_features <- input$features %||% character()
        cur_use_core <- isTRUE(input$use_core)
        cur_core     <- if (cur_use_core) (input$core_features %||% character()) else character()
        cur_mode     <- input$mode
        cur_k        <- input$k_range
        cur_panel    <- input$panel
        cur_truth   <- input$truth_column

        !setequal(cur_features, committed$features) ||
        !identical(cur_use_core, committed$use_core) ||
        !setequal(cur_core, committed$core) ||
        !identical(cur_mode, committed$mode) ||
        !identical(cur_k, committed$k_range) ||
        !identical(cur_panel, committed$panel) ||
        !identical(cur_truth, committed$truth_column)
    })

    observeEvent(TRUE, {
        committed$truth_column <- input$truth_column %||% "lymphgen"

        committed$features <- input$features %||% character()
        committed$use_core <- isTRUE(input$use_core)
        committed$core     <- if (committed$use_core) (input$core_features %||% character()) else character()
        committed$mode     <- input$mode %||% "Lenient"
        committed$k_range  <- input$k_range %||% c(k_low, k_high)
        committed$panel    <- input$panel %||% default_panel

        init_done(TRUE)
        shinyjs::enable("predict")  # allow Run Classifier on fresh load
    }, once = TRUE)


    # Build a one-row test df for a custom sample, aligned to the current model's feature set
    make_custom_test_df <- function(sample_id, selected_features, present_value, train_cols) {
        # start with zeros for all columns in the model
        row_vals <- setNames(as.list(rep(0, length(train_cols))), train_cols)

        # set selected features to present_value (only those that exist in current train set)
        selected_features <- intersect(selected_features, train_cols)
        for (f in selected_features) row_vals[[f]] <- present_value

        # return as data.frame with a rowname/sample_id column preserved later by predict code
        out <- as_tibble(row_vals)
        # ensure numeric
        out[] <- lapply(out, function(x) as.numeric(x))
        sane_sample_id <- sanitize_sample_id(sample_id, existing_ids = rownames(full_status))
        rownames(out) <- sane_sample_id
        out
    }
    sanitize_sample_id <- function(x, existing_ids = character()) {
        y <- trimws(x)
        y <- gsub("[^A-Za-z0-9._-]", "_", y)      # allow letters, numbers, . _ -
        if (!grepl("^[A-Za-z]", y)) y <- paste0("S_", y)  # ensure starts with a letter
        # avoid leading '.' (hidden names)
        if (grepl("^\\.", y)) y <- sub("^\\.", "S_", y)
        # enforce uniqueness
        if (y %in% existing_ids) {
        y <- make.unique(c(existing_ids, y), sep = "_")
        y <- tail(y, 1)
        }
        y
    }
    observe({
        if (isTRUE(dirty())) {
        shinyjs::disable("predict")
        } else {
        shinyjs::enable("predict")
        }
    })

    observe({
        session$registerDataObj("download_predictions", list(run_id = NULL), function(data, req) {
        # Extract run_id from query string
        run_id <- req$QUERY_STRING %>%
            strsplit("=", fixed = TRUE) %>%
            unlist() %>%
            .[2]

        # isolate() avoids requiring a reactive context
        stored_predictions <- isolate(prediction_store())
        preds <- stored_predictions[[run_id]]

        if (!is.null(preds)) {
            tmp <- tempfile(fileext = ".tsv")
            write_tsv(preds, file = tmp)
            return(function(res) {
            res$setHeader("Content-Type", "text/tsv")
            res$setHeader("Content-Disposition", paste0("attachment; filename=DLBCLone_predictions_", run_id, ".tsv"))
            res$sendFile(tmp)
            })
        } else {
            return(function(res) {
            res$setHeader("Content-Type", "text/plain")
            res$write(paste("Error: No prediction data found for run_id:", run_id))
            res$finish()
            })
        }
        })
    })

    observe({
        # Only run once
        if (is.null(run_log())) {
        feature_str <- paste(colnames(full_status), collapse = ",")
        default_link <- paste0(
            "?panel=", default_panel,
            "&k_min=", k_low,
            "&k_max=", k_high,
            "&features=", URLencode(feature_str)
        )
        run_id <- paste0(ncol(full_status), "_features_", as.integer(Sys.time()))

        current_store <- prediction_store()
        current_store[[run_id]] <- default_knn$predictions
        prediction_store(current_store)


        # Save predictions to disk
        # csv_path <- file.path(pred_dir, paste0("DLBCLone_predictions_", run_id, ".csv"))
        tsv_path <- file.path(pred_dir, paste0("DLBCLone_predictions_", run_id, ".tsv"))
        write_tsv(default_knn$predictions, file = tsv_path)

        run_log(data.frame(
            timestamp = format(Sys.time(), "%b %d, %Y %I:%M %p"),
            run_id = run_id,
            panel = default_panel,
            Unclass = filter(default_knn$predictions, DLBCLone_ko == "Other") %>% nrow(),
            k_range = paste(k_low, k_high, sep = " - "),
            k_used = default_knn$DLBCLone_k_best_k,
            accuracy = round(default_knn$DLBCLone_k_accuracy, 3),
            purity_threshold = round(default_knn$DLBCLone_k_purity_threshold, 3),
            n_features = ncol(full_status),
            download = paste0(
            "<a class='btn btn-sm btn-primary' href='predictions/DLBCLone_predictions_", run_id,
            ".tsv' download>Download</a>"
            ),
            Revisit = paste0("<a href='", default_link, "' target='_blank'>Link</a>"),
            stringsAsFactors = FALSE
        ))
        }
    })
    #initialize
    observe({
        if(is.null(sample_log())){
        sample_log(data.frame(
            sample_id = character(),
            prediction = character(),
            stringsAsFactors = FALSE
        ))
        }
    })

    observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!length(query)) return()

    restoring_from_query(TRUE)

    if (!is.null(query$use_core)) {
        updateCheckboxInput(session, "use_core",
        value = tolower(query$use_core) %in% c("1","true","t","yes","y"))
    }

    if (!is.null(query$panel))   updateSelectInput(session, "panel", selected = query$panel)
    if (!is.null(query$mode))    updateRadioButtons(session, "mode", selected = query$mode)
    if (!is.null(query$k_min) && !is.null(query$k_max))
        updateSliderInput(session, "k_range", value = c(as.numeric(query$k_min), as.numeric(query$k_max)))

    if (!is.null(query$features)) {
        feature_list <- strsplit(query$features, ",")[[1]]
        updateCheckboxGroupInput(session, "features", selected = feature_list)
    }

    if (!is.null(query$core)) {
        core_list <- intersect(strsplit(query$core, ",")[[1]], sort(colnames(full_status)))
        pending_core(core_list)  # <- just store; we'll apply when inputs are ready
    }

    # leave restoring_from_query(TRUE) until features are processed (see next section)
    })


    #default features (all)
    observeEvent(input$panel, {
        allowed_all <- sort(colnames(full_status))
        updateCheckboxGroupInput(
        session, "features",
        choices  = allowed_all,
        selected = if (isTRUE(restoring_from_query())) isolate(input$features) else panels[[input$panel]],
        inline   = TRUE
        )
    }, ignoreInit = TRUE)


    ## Testing permalink
    observeEvent(input$share_link, {
        # snapshot the last-trained settings
        params <- list(
        panel    = committed$panel,
        k_min    = committed$k_range[1],
        k_max    = committed$k_range[2],
        mode     = committed$mode,
        use_core = as.integer(committed$use_core),            # 0/1
        features = paste(committed$features, collapse = ",")
        )
        if (isTRUE(committed$use_core) && length(committed$core)) {
        params$core <- paste(committed$core, collapse = ",")
        }

        url <- paste0(base_url(), "?", build_query(params))

        showModal(modalDialog(
    title = "Return link to this model setup",
    easyClose = TRUE,
    footer = tagList(
        modalButton("Close"),
        tags$button(
        id = "copy_permalink_btn",
        type = "button",
        class = "btn btn-primary",
        "Copy link"
        )
    ),
    tagList(
        p("Open this URL to restore these settings in the app:"),
        tags$textarea(
        id = "permalink_text",
        style = "width:100%;height:80px;",
        url
        ),
        tags$hr(),
        p(tags$em("Tip: share this link or bookmark it to recreate the model later."))
    )
    ))

    })


    observeEvent(input$predict, {
        current_DLBCLone_KNN <- dlbclone_result()

        # Columns the model expects (training columns)
        train_df <- full_status %>%
        select(all_of(colnames(current_DLBCLone_KNN$features_df)))
        #core_feats = core_features()
        if (identical(input$predict_source, "Demo sample")) {
        # -------- DEMO path (unchanged logic) --------
        req(input$predict_sample)
        sample_features <- prototypes[input$predict_sample, , drop = FALSE] %>%
            select(all_of(colnames(current_DLBCLone_KNN$features_df)))
        #print("FEATS:")
        #print(colnames(sample_features))
        predicted <- DLBCLone_KNN_predict(
            train_df = train_df,
            test_df = sample_features,
            metadata = dlbcl_meta_clean,
            #core_features = core_feats,
            DLBCLone_KNN_out = current_DLBCLone_KNN,
            truth_column = current_truth_column()
        )

        dlbclone_pred_result(predicted)
        current_predict_id(input$predict_sample)

        # Columns present in prediction
        pred_cols <- colnames(predicted$unlabeled_predictions)

        # Class score columns are named per class; Plus the backend keeps 'Other_score' stable
        class_score_cols <- intersect(pred_cols, current_classes())
        extra_cols <- intersect(
            c("Other_score", "valid_classes", "score_ratio", "neighbor_id"),
            pred_cols
        )

        sample_id_to_show <- if (identical(input$predict_source, "Demo sample"))
            input$predict_sample else input$custom_sample_id

        sample_class <- predicted$unlabeled_predictions %>%
            filter(sample_id == sample_id_to_show) %>%
            select(all_of(c("sample_id", "DLBCLone_ko", class_score_cols, extra_cols)))

        } else {
        # -------- CUSTOM path --------
        validate(
            need(nzchar(input$custom_sample_id), "Please enter a Sample ID."),
            need(length(input$custom_selected) > 0, "Select at least one feature for your custom sample.")
        )

        test_df <- make_custom_test_df(
            sample_id = input$custom_sample_id,
            selected_features = input$custom_selected,
            present_value = input$custom_present_value,
            train_cols = colnames(current_DLBCLone_KNN$features_df)
        )
        #core_feats = core_features()
        #print(test_df)
        # Predict
        predicted <- DLBCLone_KNN_predict(
            train_df = train_df,
            test_df = test_df,
            metadata = dlbcl_meta_clean,
            #core_features = core_feats,
            DLBCLone_KNN_out = current_DLBCLone_KNN
        )

        dlbclone_pred_result(predicted)
        current_predict_id(input$custom_sample_id)

        # Columns present in prediction
        pred_cols <- colnames(predicted$unlabeled_predictions)

        # Class score columns are named per class; Plus the backend keeps 'Other_score' stable
        class_score_cols <- intersect(pred_cols, current_classes())
        extra_cols <- intersect(
            c("Other_score", "valid_classes", "score_ratio", "neighbor_id"),
            pred_cols
        )

        sample_id_to_show <- if (identical(input$predict_source, "Demo sample"))
            input$predict_sample else input$custom_sample_id

        sample_class <- predicted$unlabeled_predictions %>%
            filter(sample_id == sample_id_to_show) %>%
            select(all_of(c("sample_id", "DLBCLone_ko", class_score_cols, extra_cols)))

        }
        # Round all numeric columns to 3 decimal places
        sample_class <- sample_class %>%
        mutate(across(where(is.numeric), ~ round(.x, 3)))
        #add more details
        sample_class$panel = input$panel
        sample_class$Mode = input$mode
        sample_class$n_feats = length(input$features)
        sample_class$core = paste0(input$core_features, collapse = ",")

        # Append to sample log
        sample_class =
        sample_class %>%
        rename(prediction=DLBCLone_ko) %>%
        relocate(panel:core, .after = prediction)
        cur_log <- sample_log()


        cur_log <- bind_rows(cur_log, sample_class)
        sample_log(cur_log)
    })


    # re-generate and evaluate the classifier
    observeEvent(input$regenerate, {
        status <- full_status %>% select(all_of(input$features))
        optimize_for_other = if_else(input$mode == "Stringent", TRUE, FALSE)
        truth_col <- input$truth_column
        other_lbl <- "Other"
        #core_feats = core_features()
        core_feats = selected_core()
        meta_classes <- dlbcl_meta_clean[[truth_col]] %>% unique() %>% setdiff(NA) %>% as.character()
        meta_classes = sort(meta_classes[meta_classes!=other_lbl])

        meta_classes = c(meta_classes,other_lbl)
        updated_result <- DLBCLone_KNN(status,
        dlbcl_meta_clean,
        core_features = core_feats,
        min_k = input$k_range[1],
        max_k = input$k_range[2],
        optimize_for_other = optimize_for_other,
        truth_column = truth_col,
        truth_classes = meta_classes,
        other_class = other_lbl
        )
        dlbclone_result(updated_result)

        feature_str <- paste(input$features, collapse = ",")
        link <- paste0(
        # session$clientData$url_hostname, session$clientData$url_pathname,
        "?panel=", input$panel,
        "&sample=", input$sample,
        "&k_min=", input$k_range[1],
        "&k_max=", input$k_range[2],
        "&features=", URLencode(feature_str)
        )
        run_id <- paste0(length(input$features), "_features_", as.integer(Sys.time()))

        tsv_path <- file.path(pred_dir, paste0("DLBCLone_predictions_", run_id, ".tsv"))
        write_tsv(updated_result$predictions, file = tsv_path)

        log_df <- run_log()
        log_df <- rbind(
        log_df,
        data.frame(
            timestamp = format(Sys.time(), "%b %d, %Y %I:%M %p"),
            run_id = run_id,
            panel = input$panel,
            Unclass = filter(updated_result$predictions, DLBCLone_ko == "Other") %>% nrow(),
            k_range = paste(input$k_range[1], input$k_range[2], sep = " - "),
            k_used = updated_result$DLBCLone_k_best_k,
            accuracy = round(updated_result$DLBCLone_k_accuracy, 3),
            purity_threshold = round(updated_result$DLBCLone_k_purity_threshold, 3),
            n_features = length(input$features),
            download = paste0(
            "<a class='btn btn-sm btn-primary' href='predictions/DLBCLone_predictions_", run_id,
            ".tsv' download>Download</a>"
            ),
            Revisit = paste0("<a href='", link, "' target='_blank'>Link</a>"),
            stringsAsFactors = FALSE
        )
        )
        run_log(log_df)
        # snapshot the settings that produced this model
        committed$features <- input$features %||% character()
        committed$use_core <- isTRUE(input$use_core)
        committed$core     <- if (committed$use_core) (input$core_features %||% character()) else character()
        committed$mode     <- input$mode
        committed$k_range  <- input$k_range
        committed$panel    <- input$panel
        committed$truth_column <- input$truth_column

        # (optional) ensure Predict is enabled immediately after training
        shinyjs::enable("predict")

    })
    observeEvent(list(input$use_core, input$core_features), {
        if (!isTRUE(init_done())) return()

        cur_use_core <- isTRUE(input$use_core)
        cur_core     <- if (cur_use_core) (input$core_features %||% character()) else character()

        if (!identical(cur_use_core, committed$use_core) ||
            !setequal(cur_core, committed$core)) {
        shinyjs::disable("predict")
        }
    }, ignoreInit = TRUE)

    observeEvent(input$features, {
    allowed <- sort(input$features %||% character())

    if (isTRUE(restoring_from_query())) {
        # If we have URL-provided cores, apply them now (intersect with allowed)
        cores <- pending_core()
        if (!is.null(cores)) {
        updateSelectizeInput(session, "core_features",
                            choices  = allowed,
                            selected = intersect(cores, allowed),
                            server   = TRUE)
        pending_core(NULL)
        } else {
        # No specific cores to set; just refresh choices without clobbering selection
        updateSelectizeInput(session, "core_features",
                            choices = allowed,
                            server  = TRUE)
        }

        # We’re done restoring once features (and possibly cores) are applied
        restoring_from_query(FALSE)
        return()
    }

    # Normal path (not restoring): keep current core selection but drop any invalids
    current_core <- isolate(input$core_features) %||% character()
    new_core <- intersect(current_core, allowed)
    updateSelectizeInput(session, "core_features",
                        choices  = allowed,
                        selected = new_core,
                        server   = TRUE)
    }, ignoreInit = TRUE)


    observeEvent(input$visualize, {
        status <- full_status %>% select(all_of(input$features))
        # core = c("MYD88HOTSPOT","NOTCH1","SGK1","EZH2","NOTCH2","BCL6_SV","TET2")
        # core_features = core[core %in% input$features]
        #fc = core_features()
        fc  = selected_core()
        updated_result <- make_and_annotate_umap(
        df = status,
        metadata = dlbcl_meta_clean,
        core_features = fc
        )


        umap_result(updated_result)
    })
    output$run_log_table <- renderDT({
        datatable(run_log(), escape = FALSE, options = list(pageLength = 10))
    })
    output$sample_log_table <- renderDT({
        datatable(sample_log(), escape = FALSE, options = list(pageLength = 10))
    })

    output$heatmap <- renderPlot({
        DLBCLone_KNN_out <- dlbclone_pred_result()
        sid <- current_predict_id()

        if (is.null(DLBCLone_KNN_out) ||
            is.null(DLBCLone_KNN_out$predictions) ||
            is.null(sid) ||
            !(sid %in% DLBCLone_KNN_out$unlabeled_predictions$sample_id)) {
        plot.new()
        text(0.5, 0.5, "Run the classifier (Demo or Custom) to view neighbors.", cex = 1.5)
        return()
        }
        if(is.null(DLBCLone_KNN_out$unlabeled_neighbors)) {
        plot.new()
        text(0.5, 0.5, "No neighbors found.", cex = 1.5)
        return()
        }
        print(DLBCLone_KNN_out$unlabeled_predictions)
        print(DLBCLone_KNN_out$unlabeled_neighbors)
        nearest_neighbor_heatmap(sid, DLBCLone_KNN_out, truth_column = current_truth_column())
    })
    output$DLBCLone_KNN_plot_truth <- renderPlotly({
        result <- dlbclone_result()
        if (is.null(result)) {
        result <- default_knn
        }
        validate(need(!is.null(result), "Waiting for model output..."))
        #ggplotly(result$plot_truth, tooltip = c("sample_id", "lymphgen", "DLBCLone_ko")) %>%
        ggplotly(result$plot_truth, tooltip = c("sample_id", current_truth_column(), "DLBCLone_ko")) %>%
        plotly::layout(width = 800, height = 600)
    })
    output$DLBCLone_KNN_plot_prediction <- renderPlotly({
        result <- dlbclone_result()
        # if(is.null(result)){
        #  result = default_knn
        # }
        validate(need(!is.null(result), "Waiting for model output..."))
        #ggplotly(result$plot_predicted, tooltip = c("sample_id", "lymphgen", "DLBCLone_ko")) %>%
        ggplotly(result$plot_predicted, tooltip = c("sample_id", current_truth_column(), "DLBCLone_ko")) %>%
        plotly::layout(width = 800, height = 600)
    })

    output$predictions <- renderDT({
        result <- dlbclone_result()
        # if(is.null(result)){
        #  result = default_knn
        # }
        datatable(
        result$predictions %>%
            #mutate(
            #  top_group_score = round(top_group_score, 3),
            #  score_ratio = round(score_ratio, 5)
            #) %>%
            #select(sample_id, top_group_score, score_ratio, lymphgen, starts_with("DLBCLone")),
            select(all_of(c("sample_id", "top_group_score", "score_ratio",
                    current_truth_column())), starts_with("DLBCLone")),
        options = list(pageLength = 50, searching = TRUE)
        )
    })

    output$umap_scatterplot <- renderPlot({
        umap_out <- umap_result()
        # if(is.null(umap_out)){
        #  umap_out = default_umap
        # }
        validate(need(!is.null(umap_out), "Waiting for model output..."))
        tc = current_truth_column()
        make_umap_scatterplot(umap_out$df,  colour_by = tc)
    })

    output$alluvial <- renderPlot({
        result <- dlbclone_result()
        if (is.null(result)) {
        result <- default_knn
        }
        validate(need(!is.null(result), "Waiting for model output..."))
        make_alluvial(result,
        # count_excluded_as_other = T,
        pred_column = "DLBCLone_ko",

        label_size = 4,
        title = paste0("Optimal K=", result$DLBCLone_k_best_k, ",")
        )
    })
    }

    # Run the application
    shinyApp(ui = ui, server = server)
}
