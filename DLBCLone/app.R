#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(GAMBLR)
library(GAMBLR.predict)
lyseq_genes = read.table(paste0(here::here(),"/schmitz_lstgene_list.tsv")) %>% pull(V1)
lyseq_status = read_tsv(file=paste0(here::here(),"/lyseq_status.tsv")) %>% column_to_rownames("sample_id")
dlbcl_meta_clean = read_tsv(file=paste0(here::here(),"/dlbcl_meta_clean.tsv"))

## make demo data


prototypes = lyseq_status[c(1:12),]
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


#DLBCLone_out = DLBCLone_KNN(lyseq_status,dlbcl_meta_clean,min_k=10,max_k=10)


ui <- fluidPage(

  # Application title
  titlePanel("DLBCLone"),

  sidebarLayout(
    sidebarPanel(
      actionButton("submit_button", "Re-run", icon("refresh")),  # ✅ changed from submitButton
      helpText("When you click the button a new model will be generated using the selected features"),
      checkboxGroupInput("features", label = "Features",
                         choices = colnames(lyseq_status),
                         selected = colnames(lyseq_status))
    ),

    mainPanel(
      plotOutput("DLBCLone_KNN_plot"),
      plotOutput("alluvial")
    )
  )
)

server <- function(input, output) {

  dlbclone_result <- reactiveVal()

  # Initial run
  observe({
    initial_result <- DLBCLone_KNN(lyseq_status, dlbcl_meta_clean, plot_samples = "00-17960_CLC01670")
    dlbclone_result(initial_result)
  })

  # Run again when user clicks "Re-run"
  observeEvent(input$submit_button, {
    status = lyseq_status %>% select(all_of(input$features))
    updated_result <- DLBCLone_KNN(status, dlbcl_meta_clean, plot_samples = "00-17960_CLC01670")
    dlbclone_result(updated_result)
  })

  output$DLBCLone_KNN_plot <- renderPlot({
    result = dlbclone_result()
    validate(need(!is.null(result), "Waiting for model output..."))
    result$plot
  })

  output$alluvial <- renderPlot({
    result = dlbclone_result()
    validate(need(!is.null(result), "Waiting for model output..."))
    make_alluvial(result, pred_column = "DLBCLone_ko")
  })
}

# Run the application
shinyApp(ui = ui, server = server)
