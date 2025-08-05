#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(DT)
library(shinybusy)
library(GAMBLR)
library(ggalluvial)
library(caret)
library(GAMBLR.predict)
lyseq_genes = read.table(paste0(here::here(),"/schmitz_lstgene_list.tsv")) %>% pull(V1)
lyseq_status = read_tsv(file=paste0(here::here(),"/lyseq_status.tsv")) %>% column_to_rownames("sample_id")
full_status = read_tsv(file=paste0(here::here(),"/all_full_status.tsv")) %>% column_to_rownames("sample_id")
dlbcl_meta_clean = read_tsv(file=paste0(here::here(),"/dlbcl_meta_clean.tsv"))
dlbcl_meta_clean$cohort = "GAMBL"
# This code needs to be fixed to work anywhere
lacy_df = readxl::read_excel("~/git/Lyntegrate/data/bloodbld2019003535-suppl2.xlsx",sheet=1) %>%
  mutate(Gene= gsub("_.+","",Gene))
all_lacy_genes = c(unique(lacy_df$Gene),"MYD88HOTSPOT")
all_lacy_genes = all_lacy_genes[all_lacy_genes %in% colnames(full_status)]


lacy_df =lacy_df %>% filter(`Included in statistical analysis`=="Yes")

lacy_genes = lacy_df %>% pull(Gene) %>% unique()
lacy_genes = c(lacy_genes,"MYD88HOTSPOT")
lacy_genes = lacy_genes[lacy_genes %in% colnames(full_status)]
## make demo data



prototypes = full_status[c(1:12),]
rownames(prototypes) = c("BN2_1","BN2_2","ST2_1","ST2_2","EZB_1","EZB_2","MCD_1","MCD_2","N1_1","N1_2","Other_1","Other_2")
prototypes[]=0
prototypes["BN2_1",c("BCL6_SV","TNFAIP3","KLF2")]=1
prototypes["BN2_2",c("NOTCH1","BCL6","SPEN","UBE2A")]=1
prototypes["ST2_1",c("SGK1","ZFP36L1","ACTG1","STAT3","CD83")]=1
prototypes["ST2_2",c("TET2","ITPKB","NFKBIA","DUSP2","JUNB")]=1
prototypes["MCD_1",c("MYD88HOTSPOT","PIM1","ETV6")]=1
prototypes["MCD_2",c("CD79B","PIM1","BTG1","GRHPR")]=1
prototypes["EZB_1",c("CREBBP","BCL2_SV","KMT2D")]=1
prototypes["EZB_2",c("EZH2","TNFRSF14","IRF8")]=1
prototypes["N1_1",c("NOTCH1","ATM")] = 1
prototypes["N1_2",c("NOTCH1","ID3")] = 1
prototypes["Other_1",c("CD83","TP53","BTK","NFKBIZ","STAT6")]=1
prototypes["Other_2",c("RRAGC","P2RY8","CDKN2A","MS4A1","B2M")]=1

full_status = bind_rows(full_status,prototypes)


default_umap  <- make_and_annotate_umap(df = full_status,
                                        metadata = dlbcl_meta_clean)
k_low = 5
k_high = 20
default_knn <- DLBCLone_KNN(full_status, dlbcl_meta_clean,
                            min_k = k_low,
                            max_k = k_high,
                            #core_features = core_features,
                            predict_unlabeled = T)
impact_genes = c("MYD88HOTSPOT","ABL1","BCL2","CEBPA","ETV6","HGF","JUN","MSH2","PHF6","RPTOR","SRSF2","ZRSR2","ACTG1","BCL6","CHEK1","EZH2","HIF1A","KDM5A","MSH6","PIGA","RRAGC","STAG1","H1B","AKT1","BCOR","CHEK2","FAM46C","HIST1H1B","KDM5C","MTOR","PIK3C2G","RTEL1","STAG2","H1-2","AKT2","BCORL1","CIC","FANCA","HIST1H1C","KDM6A","MUTYH","PIK3C3","RUNX1","STAT3","H1D","AKT3","BCR","CIITA","FANCC","HIST1H1D","KDR","MYC","PIK3CA","RUNX1T1","STAT5A","H1E","ALK","BIRC3","CRBN","FANCD2","HIST1H1E","KEAP1","MYCL1","PIK3CG","SAMHD1","STAT5B","H2AC","ALOX12B","BLM","CREBBP","FAS","HIST1H2AC","KIT","MYCN","PIK3R1","SDHA","STAT6","H2AG","AMER1","BRAF","CRKL","FAT1","HIST1H2AG","KMT2A","MYD88","PIK3R2","SDHB","STK11","H2AL","APC","BRCA1","CRLF2","FBXO11","HIST1H2AL","KMT2B","NBN","PIM1","SDHC","SUFU","H2AM","AR","BRCA2","CSF1R","FBXW7","HIST1H2AM","KMT2C","NCOR1","PLCG1","SDHD","SUZ12","H2BC","ARAF","BRD4","CSF3R","FGF19","HIST1H2BC","KMT2D","NCOR2","PLCG2","SETBP1","SYK","H2BC5","ARHGEF28","BRIP1","CTCF","FGF3","HIST1H2BD","KRAS","NCSTN","PMS2","SETD1A","TBL1XR1","H2BG","ARID1A","BTG1","CTNNB1","FGF4","HIST1H2BG","KSR2","NF1","PNRC1","SETD1B","TBX3","H2BJ","ARID1B","BTK","CUX1","FGFR1","HIST1H2BJ","LCK","NF2","POT1","SETD2","TERT","H2BK","ARID2","CALR","CXCR4","FGFR2","HIST1H2BK","LMO1","NFE2","PPP2R1A","SETD3","TET1","H2BO","ARID3A","CARD11","CYLD","FGFR3","HIST1H2BO","LTB","NFE2L2","PRDM1","SETD4","TET2","H3C2","ARID3B","CASP8","DAXX","FGFR4","HIST1H3B","MALT1","NKX2-1","PRKAR1A","SETD5","TET3","H3C8","ARID3C","CBFB","DDR2","FLCN","HIST1H3G","MAP2K1","NOTCH1","PTCH1","SETD6","TGFBR2","ARID4A","CBL","DDX3X","FLT1","HLA-A","MAP2K2","NOTCH2","PTEN","SETD7","TNFAIP3","ARID4B","CCND1","DIS3","FLT3","HNF1A","MAP2K4","NOTCH3","PTPN1","SETD8","TNFRSF14","ARID5A","CCND2","DNMT3A","FLT4","HRAS","MAP3K1","NOTCH4","PTPN11","SETDB1","TOP1","ARID5B","CCND3","DOT1L","FOXL2","ID3","MAP3K13","NPM1","PTPN2","SETDB2","TP53","ASXL1","CCNE1","DTX1","FOXO1","IDH1","MAP3K14","NRAS","RAD21","SF3B1","TP63","ASXL2","CD274","DUSP22","FOXP1","IDH2","MAPK1","NSD1","RAD50","SGK1","TRAF2","ATM","CD28","EED","FURIN","IGF1","MAPK3","NT5C2","RAD51","SH2B3","TRAF3","ATP6AP1","CD58","EGFR","FYN","IGF1R","MCL1","NTRK1","RAD51B","SMAD2","TRAF5","ATP6V1B2","CD79A","EGR1","GATA1","IGF2","MDM2","NTRK2","RAD51C","SMAD4","TSC1","ATR","CD79B","EP300","GATA2","IKBKE","MDM4","NTRK3","RAD51D","SMARCA4","TSC2","ATRX","CDC73","EP400","GATA3","IKZF1","MED12","P2RY8","RAD52","SMARCB1","TSHR","ATXN2","CDH1","EPHA3","GNA11","IKZF3","MEF2B","PAK7","RAD54L","SMARCD1","TYK2","AURKA","CDK12","EPHA5","GNA12","IL7R","MEN1","PALB2","RAF1","SMC1A","U2AF1","AURKB","CDK4","EPHA7","GNA13","INPP4B","MET","PARP1","RARA","SMC3","U2AF2","AXIN1","CDK6","EPHB1","GNAQ","IRF1","MGA","PAX5","RB1","SMG1","UBR5","AXL","CDK8","ERBB2","GNAS","IRF4","MGAM","PBRM1","REL","SMO","VAV1","B2M","CDKN1B","ERBB3","GNB1","IRF8","MITF","PCBP1","RET","SOCS1","VAV2","BACH2","CDKN2A","ERBB4","GRIN2A","IRS2","MLH1","PDCD1","RHOA","SOX2","VHL","BAP1","CDKN2Ap14ARF","ERG","GSK3B","JAK1","MOB3B","PDGFRA","RICTOR","SP140","WHSC1","BARD1","CDKN2Ap16INK 4A","ESCO2","HDAC1","JAK2","MPEG1","PDGFRB","RNF43","SPEN","WT1","BCL10","CDKN2B","ESR1","HDAC4","JAK3","MPL","PDPK1","ROBO1","SPOP","XBP1","BCL11B","CDKN2C","ETNK1","HDAC7","JARID2","MRE11A","PDS5B","ROS1","SRC","XPO1")
lymphgen_genes = c( "ACTB"  ,     "ACTG1"   ,   "BCL10"   ,   "BCL2"   ,
                   "BCL2L1"  ,   "BCL6"     ,  "BTG1"  ,     "BTG2"    ,
                    "CD70"   ,    "CD79B"  ,    "CD83"   ,    "CREBBP"    ,
                   "DDX3X"  ,    "DTX1"     ,  "DUSP2"    ,  "EDRF1"     ,
                   "EIF4A2"    , "EP300"  ,    "ETS1"   ,    "ETV6"  ,
                    "EZH2",      "FOXC1"   ,   "GRHPR"  ,    "HIST1H1B"  ,
                    "HIST1H2BC" , "HLA-A"  ,    "HLA-B"  ,    "ID3"    ,
                    "IRF8"   ,    "ITPKB" ,     "JUNB"       ,"KLF2"  ,
                    "KLHL14"    , "KMT2D"  ,    "MBTPS1" ,    "MED16"   ,
                    "MEF2B" ,     "MPEG1"  ,    "MYD88HOTSPOT", "NFKBIA"    ,
                    "NOL9" ,      "NOTCH1" ,    "NOTCH2"    , "OSBPL10"   ,
                    "PIM1"  ,     "PIM2"  ,     "PRDM1"  ,    "PRRC2A",
                    "PPRC2C" ,    "RFTN1" ,     "SEC24C"  ,   "SETD1B" ,
                    "SGK1" ,      "SOCS1" ,     "SPEN"  ,     "STAT3"  ,
                   "TBL1XR1",    "TET2"  ,     "TNFAIP3"  ,  "TNFRSF14" ,
                   "UBE2A" ,     "WEE1"   ,    "ZFP36L1"  , "BCL2_SV", "BCL6_SV" )
lyseq_genes = sort(colnames(lyseq_status)[colnames(lyseq_status) %in% colnames(full_status)])
panels = list(everything=sort(colnames(full_status)),
              Lacy_full = all_lacy_genes,
              Lacy_selected = lacy_genes,
              Lymphgen_not_CNV = lymphgen_genes[lymphgen_genes %in% colnames(full_status)],
              Lymphplex = c("BCL2",	"BCL6",	"ARID1A",	"B2M","BTG1",	"BTG2",	"CCND3",
                            "CD70",	"CD79B",	"CIITA",	"CREBBP",	"DDX3X",	"DTX1",
                            "DUSP2",	"EP300",	"EZH2",	"FAS",	"GNA13",	"IRF4",
                            "IRF8",	"KMT2D",	"MPEG1",	"MYD88",	"NOTCH1",	"NOTCH2",	"PIM1","PRDM1","SGK1",	"SOCS1",
                            "STAT3",	"STAT6"	,"TBL1XR1",	"TET2","TNFAIP3","TNFRSF14",	"ZFP36L1"),
              Lyseq = lyseq_genes,
              Lyseq_no_TP53 = lyseq_genes[lyseq_genes != "TP53"],
              MSK_impact =impact_genes[impact_genes %in% colnames(full_status)],
              mini_set1 = c("KLHL14","TMEM30A","RHOA","HIST1H1D","CARD11","CREBBP","FOXO1",
                            "SPEN","CD79B","SETD1B","DTX1","STAT3","HIST1H1B","TET2",
                            "S1PR2","EZH2","STAT6","EP300","ACTG1","BCL2",
                            "ETS1","PIM1","NOTCH1","SGK1","NOTCH2"),
              mini_set2 = c( "BCL2","TET2", "ZFP36L1","EP300","HLA-A" ,    "BCL10"  ,   "HIST1H2BC",
                             "CD79B"  ,   "IRF8"    ,  "DTX1"   ,   "HIST1H1B",  "ETS1"  ,    "CREBBP" ,   "CARD11" ,
                             "CXCR4"  ,   "TNFAIP3",   "BTG1"  ,    "KLHL14"   , "KMT2D" ,    "STAT6"  ,   "NOTCH1" ,
                             "B2M"  ,     "NFKBIA"  ,  "ATM"   ,    "TBL1XR1" ,  "P2RY8" ,    "PIM1" ,     "SGK1" ,
                             "NOTCH2"))

demo_samples = rownames(prototypes)

ui <- fluidPage(
  add_busy_spinner(spin = "fading-circle", color = "#003366",position = "full-page"),
  tags$div(id = "loading-text", style = "color:#003366; font-weight:bold; display:none;",
           "Running model... please wait."),
  # Application title
  titlePanel("DLBCLone"),

  sidebarLayout(
    sidebarPanel(

      actionButton("visualize", "Update UMAP", icon("refresh")),
      actionButton("regenerate", "Update Classifier", icon("refresh")),
      helpText("Visualize selected features or create a custom classifier using them"),
      selectInput("panel",
                  label="Pre-defined gene panel",
                  choices = names(panels),
                  selected = "everything"),
      #selectInput("sample",
      #            label="Label a simulated sample",
      #            choices = demo_samples),
      checkboxGroupInput("features",
                         label = "Features",inline = T,
                         choices = sort(colnames(full_status)),
                         selected = sort(colnames(full_status))),
      sliderInput("k_range",label = "Range of K values to try",min = 5, max = 50,step = 5,value = c(k_low,k_high))

    ),

    mainPanel(
      tabsetPanel(
        tabPanel("UMAP", plotOutput("umap_scatterplot")),
        tabPanel("UMAP (Lymphgen)", plotlyOutput("DLBCLone_KNN_plot_truth")),
        tabPanel("UMAP (DLBCLone)", plotlyOutput("DLBCLone_KNN_plot_prediction")),
        tabPanel("Results overview", plotOutput("alluvial",height = "800px", width="700px")),
        tabPanel("All DLBCLone assignments",downloadButton("downloadData", "Download"), DTOutput("predictions"))
      )
    )
  )
)

server <- function(input, output, session) {

  dlbclone_result <- reactiveVal(isolate(default_knn))
  umap_result <- reactiveVal(isolate(default_umap))

  observeEvent(input$panel, {
    updateCheckboxGroupInput(
      session,
      inputId = "features",
      choices = sort(colnames(full_status)),
      selected = panels[[input$panel]], inline=T
    )
    #umap_result()$df
  })


  # Run again when user clicks "Re-run"

  observeEvent(input$regenerate, {
    status = full_status %>% select(all_of(input$features))
    #core = c("MYD88HOTSPOT","NOTCH1","SGK1","EZH2","NOTCH2","BCL6_SV","TET2")
    #core_features = core[core %in% input$features]
    updated_result <- DLBCLone_KNN(status, dlbcl_meta_clean,
                                   min_k = input$k_range[1],
                                   max_k = input$k_range[2],
                                   #core_features = core_features,
                                   predict_unlabeled = T)
    dlbclone_result(updated_result)
  })

  observeEvent(input$visualize, {
    status = full_status %>% select(all_of(input$features))
    #core = c("MYD88HOTSPOT","NOTCH1","SGK1","EZH2","NOTCH2","BCL6_SV","TET2")
    #core_features = core[core %in% input$features]
    updated_result <- make_and_annotate_umap(df = status,
                                             metadata = dlbcl_meta_clean)


    umap_result(updated_result)
  })

  output$DLBCLone_KNN_plot_truth <- renderPlotly({
    result = dlbclone_result()
    if(is.null(result)){
      result = default_knn
    }
    validate(need(!is.null(result), "Waiting for model output..."))
    ggplotly(result$plot_truth,tooltip = c("sample_id", "lymphgen", "DLBCLone_ko"))  %>%
      layout(width = 1000, height = 600)
  })
  output$DLBCLone_KNN_plot_prediction <- renderPlotly({
    result = dlbclone_result()
    #if(is.null(result)){
    #  result = default_knn
    #}
    validate(need(!is.null(result), "Waiting for model output..."))
    ggplotly(result$plot_predicted,tooltip = c("sample_id", "lymphgen", "DLBCLone_ko")) %>%
      layout(width = 1000, height = 600)
  })

  output$predictions <- renderDT({
    result = dlbclone_result()
    #if(is.null(result)){
    #  result = default_knn
    #}
    datatable(result$predictions %>%
                mutate(top_group_score=round(top_group_score,3),
                       score_ratio=round(score_ratio,5)) %>%
                select(sample_id,top_group_score,score_ratio,lymphgen,starts_with("DLBCLone"))
                ,options=list(pageLength=50,searching=TRUE))
  })

  output$umap_scatterplot <- renderPlot({

    umap_out = umap_result()
    #if(is.null(umap_out)){
    #  umap_out = default_umap
    #}
    validate(need(!is.null(umap_out), "Waiting for model output..."))
    make_umap_scatterplot(umap_out$df)

  })

  output$alluvial <- renderPlot({
    result = dlbclone_result()
    if(is.null(result)){
      result = default_knn
    }
    validate(need(!is.null(result), "Waiting for model output..."))
    make_alluvial(result,
                  #count_excluded_as_other = T,
                  pred_column = "DLBCLone_ko",
                  label_size=4,
                  title=paste0("Optimal K=",result$DLBCLone_k_best_k, ","))
  })
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("DLBCLone_predictions_",length(input$features), "_features.csv", sep = "")
    },
    content = function(file) {
      write.csv(dlbclone_result()$predictions, file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
