#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library(tidyverse)
library(ComplexHeatmap)
library(shiny)
library(shinyjs)
library(DT)
library(ggrepel)
library(shinybusy)
library(GAMBLR)
library(ggalluvial)
library(plotly)
library(caret)
library(GAMBLR.predict)

lyseq_genes <- c(
  "TNFRSF14", "NOL9", "SPEN", "ID3", "ARID1A", "RRAGC",
  "BCL10", "CD58", "NOTCH2", "BTG2", "ITPKB", "XPO1",
  "DUSP2", "CXCR4", "SF3B1", "OSBPL10", "MYD88", "RHOA",
  "NFKBIZ", "TBL1XR1", "KLHL6", "BCL6", "TET2", "IRF4",
  "CD83", "HIST1H2BC", "HIST1H1D", "HIST1H1B", "HLA-A",
  "HLA-B", "PIM1", "CCND3", "NFKBIE", "TMEM30A", "PRDM1",
  "SGK1", "TNFAIP3", "CARD11", "ACTB", "BRAF", "EZH2",
  "TOX", "MYC", "CDKN2A", "GRHPR", "NOTCH1", "FAS", "WEE1",
  "MPEG1", "MS4A1", "ATM", "ETS1", "ETV6", "KMT2D", "STAT6",
  "BTG1", "DTX1", "SETD1B", "FOXO1", "NFKBIA", "ZFP36L1",
  "B2M", "MAP2K1", "CREBBP", "CIITA", "SOCS1", "IL4R", "IRF8",
  "TP53", "STAT3", "CD79B", "GNA13", "ACTG1", "KLHL14", "BCL2",
  "TCF3", "CD70", "S1PR2", "JUNB", "KLF2", "MEF2B", "EP300",
  "P2RY8", "TMSB4X", "DDX3X", "PIM2", "BTK", "UBE2A", "MYD88HOTSPOT", "BCL2_SV", "BCL6_SV"
)

full_status <- read_tsv(file = paste0(here::here(), "/all_full_status.tsv")) %>% column_to_rownames("sample_id")

dlbcl_meta_clean <- read_tsv(file = paste0(here::here(), "/dlbcl_meta_clean.tsv"))

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



prototypes <- full_status[c(1:12), ]
rownames(prototypes) <- c("BN2_1", "BN2_2", "ST2_1", "ST2_2", "EZB_1", "EZB_2", "MCD_1", "MCD_2", "N1_1", "N1_2", "Other_1", "Other_2")
prototypes[] <- 0
maxval = max(full_status)
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

#full_status <- bind_rows(full_status, prototypes)

default_panel <- "everything"
default_umap <- make_and_annotate_umap(
  df = full_status,
  metadata = dlbcl_meta_clean
)
k_low <- 5
k_high <- 20
default_knn <- DLBCLone_KNN(full_status, dlbcl_meta_clean,
  min_k = k_low,
  max_k = k_high
)
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
  Coyle = sort(colnames(full_status)),
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
  `Lyseq (no TP53)` = lyseq_genes[lyseq_genes != "TP53"],
  `MSK Impact` = impact_genes[impact_genes %in% colnames(full_status)]
)

demo_samples <- rownames(prototypes)

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
      actionButton("regenerate", "Update Classifier", icon("refresh")),
      radioButtons("mode","Classification Mode",choices=c("Stringent","Lenient"),selected="Lenient",inline=TRUE),
      helpText("Visualize selected features or create a custom classifier using them"),
      selectInput("panel",
        label = "Pre-defined gene panel",
        choices = names(panels),
        selected = "everything"
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
      #selectInput("predict_sample",
      #            label="Predict class for a simulated sample",
      #            choices = demo_samples,
      #            selected = demo_samples[1]),
      actionButton("predict", "Run Classifier", icon("refresh")),
      helpText("Use this model to predict the class of the selected sample"),
      checkboxGroupInput("features",
        label = "Features", inline = T,
        choices = sort(colnames(full_status)),
        selected = sort(colnames(full_status))
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
      sliderInput("k_range", label = "Range of K values to try", min = 5, max = 50, step = 5, value = c(k_low, k_high))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("UMAP", plotOutput("umap_scatterplot",height = "800px",width= "700px")),
        tabPanel("UMAP (Lymphgen)", plotlyOutput("DLBCLone_KNN_plot_truth")),
        tabPanel("UMAP (DLBCLone)", plotlyOutput("DLBCLone_KNN_plot_prediction")),
        tabPanel("Results overview", plotOutput("alluvial", height = "800px", width = "700px")),
        tabPanel("All DLBCLone assignments", downloadButton("downloadData", "Download"), DTOutput("predictions")),
        tabPanel("Run Log", DTOutput("run_log_table")),
        tabPanel("Sample Neighbors", plotOutput("heatmap")),
        tabPanel("Sample Log", DTOutput("sample_log_table"))
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

  # Compute the effective feature set used for training/UMAP
  core_features <- reactive({
    base <- req(input$features)
    core = NULL
    if (isTRUE(input$use_core) && length(input$core_features)) {
      core <- input$core_features
    }
    core
  })

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
    session$registerDataObj("download_predictions", list(run_id = NULL), function(data, req) {
      # Extract run_id from query string
      run_id <- req$QUERY_STRING %>%
        strsplit("=", fixed = TRUE) %>%
        unlist() %>%
        .[2]

      # ⚠️ isolate() avoids requiring a reactive context
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

    if (!is.null(query$panel)) {
      updateSelectInput(session, "panel", selected = query$panel)
    }

    if (!is.null(query$sample)) {
      updateSelectInput(session, "sample", selected = query$sample)
    }

    if (!is.null(query$k_min) && !is.null(query$k_max)) {
      updateSliderInput(session, "k_range", value = c(as.numeric(query$k_min), as.numeric(query$k_max)))
    }

    if (!is.null(query$features)) {
      feature_list <- strsplit(query$features, ",")[[1]]
      updateCheckboxGroupInput(session, "features", selected = feature_list)
    }
  })
  #default features (all)
  observeEvent(input$panel, {
    updateCheckboxGroupInput(
      session,
      inputId = "features",
      choices = sort(colnames(full_status)),
      selected = panels[[input$panel]], inline = T
    )
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
        DLBCLone_KNN_out = current_DLBCLone_KNN
      )
      
      dlbclone_pred_result(predicted)
      current_predict_id(input$predict_sample)
      #print(names(predicted))
      #print(predicted$unlabeled_predictions)
      sample_class <- predicted$unlabeled_predictions %>%
        filter(sample_id == input$predict_sample) %>%
        #select(sample_id, DLBCLone_ko, EZB:Other_score, top_group_score, score_ratio, neighbor_id)
        select(sample_id, DLBCLone_ko, EZB:Other_score, score_ratio)

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

      sample_class <- predicted$unlabeled_predictions %>%
        filter(sample_id == input$custom_sample_id) %>%
        #select(sample_id, DLBCLone_ko, EZB:Other_score, top_group_score, score_ratio, neighbor_id)
        select(sample_id, DLBCLone_ko, EZB:Other_score, score_ratio)


      
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
    core_feats = core_features()
    # core = c("MYD88HOTSPOT","NOTCH1","SGK1","EZH2","NOTCH2","BCL6_SV","TET2")
    # core_features = core[core %in% input$features]
    updated_result <- DLBCLone_KNN(status,
      dlbcl_meta_clean,
      core_features = core_feats,
      min_k = input$k_range[1],
      max_k = input$k_range[2],
      optimize_for_other = optimize_for_other
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
  })

  observeEvent(input$features, {
    # Allowed set = currently checked Features
    allowed <- sort(if (is.null(input$features)) character() else input$features)

    # Drop any core items no longer allowed
    new_core <- intersect(isolate(input$core_features), allowed)

    # Update the core picker: limit choices to allowed; keep only valid selections
    updateSelectizeInput(
      session, "core_features",
      choices  = allowed,
      selected = new_core,
      server   = TRUE
    )
  })

  observeEvent(input$visualize, {
    status <- full_status %>% select(all_of(input$features))
    # core = c("MYD88HOTSPOT","NOTCH1","SGK1","EZH2","NOTCH2","BCL6_SV","TET2")
    # core_features = core[core %in% input$features]
    fc = core_features()
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
    nearest_neighbor_heatmap(sid, DLBCLone_KNN_out)
  })
  output$DLBCLone_KNN_plot_truth <- renderPlotly({
    result <- dlbclone_result()
    if (is.null(result)) {
      result <- default_knn
    }
    validate(need(!is.null(result), "Waiting for model output..."))
    ggplotly(result$plot_truth, tooltip = c("sample_id", "lymphgen", "DLBCLone_ko")) %>%
      layout(width = 800, height = 600)
  })
  output$DLBCLone_KNN_plot_prediction <- renderPlotly({
    result <- dlbclone_result()
    # if(is.null(result)){
    #  result = default_knn
    # }
    validate(need(!is.null(result), "Waiting for model output..."))
    ggplotly(result$plot_predicted, tooltip = c("sample_id", "lymphgen", "DLBCLone_ko")) %>%
      layout(width = 800, height = 600)
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
        select(sample_id, top_group_score, score_ratio, lymphgen, starts_with("DLBCLone")),
      options = list(pageLength = 50, searching = TRUE)
    )
  })

  output$umap_scatterplot <- renderPlot({
    umap_out <- umap_result()
    # if(is.null(umap_out)){
    #  umap_out = default_umap
    # }
    validate(need(!is.null(umap_out), "Waiting for model output..."))
    make_umap_scatterplot(umap_out$df)
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
