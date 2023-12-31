#increase filesize limit
options(shiny.maxRequestSize=10000*1024^2) #max 10gb upload

#load libraries
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(purrr)
library(stringr)
library(Seurat)
library(DT)
library(shinyalert)
library(ggplot2)


#########           UI START          ##########

# fully defined ui
ui <- dashboardPage(
      # webapp layout #
      dashboardHeader(title = "RDS Visualizer v1.1", 
                        tags$li(
                            class = "dropdown", 
                            style = "padding-top: 0px;",
                            tags$a(target="_blank", href = "https://github.com/rkafrawi/RDS_Vis_v1_1/tree/main", "Github Repo", .noWS = "outside")
                      )
  ),
  menu_bar <- dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Data Summary", tabName="input_summ", icon=icon("table")),
      menuItem("Feature Plot", tabName = "input", icon = icon("chart-bar")),
      menuItem("Dim Plot", tabName = "dim", icon = icon("cubes")),
      menuItem("Violin Plot", tabName = "violin", icon = icon("record-vinyl")),
      menuItem("Help Page", tabName = "help", icon = icon("question"))
    )
  ),
  dashboardBody(
    tabItems(
      tab_home <- tabItem(tabName = "home",
                          h2("Home Page"),
                          br(),
                          br(),
                          p("Welcome to the RDS Visualizer Tool!"),
                          br(),
                          p("This web app is designed for the visualization and exploration of single-cell RNA-seq data contained in Seurat objects. 
                             It will provide various plots and features to help you analyze and gain insights from your data. 
                             This webapp has been optimized/structured a DNAnexus build, which circumvents the instance size restrictions of shinyapps.io's base plan."),
                          br(),
                          p("This web app currently consists of four pages of note: a Data Summary page, a Feature Plot page, a Dim Plot page, and a Violin Plot page. 
                            As a general rule of thumb, these pages will not visualize any data if no .rds file has been provided by the user in the data summary page.
                            Note also that the drop downs on the Dim Plot and Violin Plot pages will not display updated categories until a file has been provided by the user."),
                          br(),
                          #embed href in paragraph sentence
                          p("Click ", a("here", target="_blank", href="https://raw.githubusercontent.com/rkafrawi/RDS_Vis_v1_1/main/docs/RDS_Visualizer_user_guide.pdf", .noWS = "outside"), " for an in depth user guide.", .noWS = c("after-begin", "before-end")),
      ),
      tab_input <- tabItem(tabName = "input_summ",
                           h2("Data Summary Page"),
                           br(),
                           #file upload
                           p("Click on the browse button and select the .rds file you wish to generate visualizations for (Max 10 GB)."),
                           fileInput("seuratFile", 
                                     "Upload RDS File",
                                     accept = ".rds"),
                           br(),
                           
                           tabsetPanel(
                             tabPanel("Metadata Overview",
                                      br(),
                                      dataTableOutput("table_cat")
                                      ),
                             tabPanel("Metadata",
                                      br(),
                                      dataTableOutput("table_meta"))
                             ),
                           
                                      
      ),
      tab_inputfeatures <- tabItem(tabName = "input",
                                   h2("Feature Plot Page"),
                                   br(),
                                   p("Once an .rds file has been uploaded, the top of this page will populate with sample genes."),
                                   br(),
                                   p("Try some of these sample genes below:"),
                                   br(),
                                   verbatimTextOutput("sampleGenes"),
                                   br(),
                                   p("Enter the feature of interest and then navigate between the UMAP and PCA subtabs to generate a feature plot for each respective reduction. Note that gene names are case sensitive."),
                                   br(),
                                   textInput("featuresInput", "Enter Gene Name:"),
                                   br(),
                                   #umap for multiple visualization tabs
                                   tabsetPanel(
                                     tabPanel("UMAP",
                                              br(),
                                              actionButton("plotButton_UMAP","Generate UMAP Plot"),
                                              br(),
                                              plotOutput("featurePlotUMAP", height = 600, width = 600)
                                              
                                     ),
                                     tabPanel("PCA",
                                              br(),
                                              actionButton("plotButton_PCA","Generate PCA Plot"),
                                              br(),
                                              plotOutput("featurePlotPCA", height = 600, width = 600))
                                     
                                   )
                                   
      ),
      # tab_input,
      # tab_feature,
      tab_dim <- tabItem(tabName = "dim",
                         h2("Dim Plot"),
                         br(),
                         p("This page generates DimPlots that can be split by different groups 
                           found in a dropdown that updates based on the metadata columns in 
                           the uploaded Seurat object. By default, the Dimplot will not be split by a group. 
                           To split the visualization by a group, toggle the checkbox and select a group from the dropdown."),
                         br(),
                         checkboxInput("splitToggle", "Split Dim Plot", value = FALSE),
                         
                         br(),
                         selectInput("variableInput", label = "Select a Group:",
                                     choices = "Input File For Dropdown Options"),
                         br(),
                         actionButton("plotButton_Dim","Generate DimPlot"),
                         br(),
                         plotOutput("featurePlotDim", height = 600, width = 900)
                         
      ),
      tab_violin <- tabItem(tabName = "violin",
                            h2("Violin Plot Page"),
                            br(),
                            p("Enter the feature of interest and then select a category to split the plot by:"),
                            br(),
                            textInput("featuresInput_vln", "Enter Gene Name:"),
                            br(),
                            selectInput("violinInput", label = "Select a Group:",
                                        choices = "Input File For Dropdown Options"),
                            br(),
                            actionButton("plotButton_violin","Generate Violin Plot"),
                            br(),
                            plotOutput("violinplot", height = 600, width = 900)
      ),
      
      tab_help <- tabItem(tabName = "help",
                          h2("Help Page"),
                          br(),
                          p("Below are a list of Frequently Asked Questions (FAQs)."),
                          br(),
                          accordion(
                            id = "helpAccordion",
                            accordionItem(
                              title = "Why is my input file not loading?",
                              "Check if the file path is correct and if the file format is compatible with the application. This webapp currently only supports .rds files with an upper size limitation of 10 GB."
                            ),
                            accordionItem(
                              title = "How do I download these figures?",
                              "To download the figures, you can right-click on each figure and choose the 'Save Image As' option. Then, select your desired location on your computer to save the image."
                            ),
                            accordionItem(
                              title = "Is there any way to view tabular representations of the visualizations?",
                              "Currently, the server does not support table generation to supplement each visual. The only way to inspect tabular representations of the uploaded data is in the Data Summary page."
                            ),
                            accordionItem(
                              title = "Why are only some of my categorical variables showing in my dropdown options?",
                              "For the sake of visual parity, only categorical variables containing less than 5 levels are included from the uploaded .rds file."
                            )
                          ),
                          br(),
                          br(),
                          br(),
                          p("Click ", a("here", target="_blank", href="https://satijalab.org/seurat/articles/pbmc3k_tutorial.html", .noWS = "outside"), " to reference the clustering guide used to build this web app.", .noWS = c("after-begin", "before-end")),
                          br(),
                          p("Sample seurat objects can be found ", a("here", target="_blank",  href="https://www.dropbox.com/home/Rizky/RShiny/rds/BMS_DGKi", .noWS = "outside"),  ". Feel free to experiment with these example datasets to familiarize yourself with the web app workflow!", .noWS = c("after-begin", "before-end")),
                          
      )
    ), 
    div(class = "hr-container",
        br(),
        hr(class = "hr-line")),
    # Footer
    div(class = "foot",
        tags$footer(
          class = "footer",
          p("Developed for use of the Belfer Center Research Team by the BxMD Group at DFCI."),
          br(),
        ),
        tags$img(
          src ="https://raw.githubusercontent.com/rkafrawi/RDS_Vis_v1_1/main/docs/bxMD_logo.jpeg",
          style="display: block; margin-left: auto; margin-right: auto; width:35%; height:15%",
        )
        ),
    
    
    
    # CSS Styling
    tags$style(HTML("
      .hr-container {
        text-align: center;
      }
      
      .hr-line {
        border-top: 1px solid #ccc;
        margin: 20px 0;
      }
      .img {
    
      }
      .footer {
        padding: 10px;
        text-align: center;
      }
      .placeholder {
      text-align: center;
      }
      .full-width-cols {
      overflow:auto;
      }
    "))
  )
)

#########           UI END          ##########


#########           SERVER START          ##########


#define server
server<- function(input, output, session) {
  
  ### Define reactive values for global vars...next ver ###
  # seuratData <- reactiveVal()
  
  ### input logic ###
  
  # Read Seurat object
  seuratData <- reactive({
    inFile <- input$seuratFile 
    if (is.null(inFile))
      return(NULL)
    
    #read Seurat object from selected file
    seurat_obj <- readRDS(inFile$datapath)
    
    return(seurat_obj)
  })
  
  # Render table-categorical only
  output$table_cat <- DT::renderDataTable({
    #retrieve seurat object from reactive expression in above code block
    seurat_data <- seuratData() 
    if (is.null(seurat_data))
      return(NULL)
    
    #get metadata
    meta_cat <- as.data.frame(seurat_data@meta.data)
    
    # include these collumns if they dont exist
    columns_to_check <- c("sample", "cytotoxic", "exact_subclonotype_id")
    
    # Find and convert columns that match the specified names
    for (col in columns_to_check) {
      if (col %in% colnames(meta_cat) && !is.factor(meta_cat[[col]])) {
        meta_cat[[col]] <- as.factor(meta_cat[[col]])
      }
    }
    
    #filter metadata by categorical
    categorical_metadata <- meta_cat[, sapply(meta_cat, is.factor)]
    
    #table of cat metadata
    datatable(categorical_metadata, options=list(scrollX=TRUE))
  })
  
  # Render table-metadata
  output$table_meta <- DT::renderDataTable({
    #retrieve seurat object from reactive expression in above code block
    seurat_data <- seuratData() 
    if (is.null(seurat_data))
      return(NULL)
    # Metadata
    seurat_df <- as.data.frame(seurat_data@meta.data)
    datatable(seurat_df, options=list(scrollX=TRUE))
  })
  
  ### featureplot logic ###
  #umap featureplot logic
  observeEvent(input$plotButton_UMAP, {
    req(seuratData(),input$featuresInput) #wait for file input and feature input
    seurat_data <- seuratData()
    if (is.null(seurat_data))
      output$errormessage <- renderText("Null file")
    
    # Split the input string @ commas into individual gene names using regex
    #[[1]] extracts the first element in list of gene names, which gives us a chr vector 
    #containing the individaul gene names as separate elements
    gene_names <- strsplit(input$featuresInput, ",\\s*")[[1]] 
    gene_names <- trimws(gene_names)  # Trim leading and trailing spaces from gene names
    
    # Check if any of the requested genes are missing
    missing_genes <- setdiff(gene_names, rownames(seurat_data@assays$RNA@data))
    
    if (length(missing_genes) > 0) {
      # output$nofeaturefound <- renderText(paste("The following genes were not found:", paste(missing_genes, collapse = ", "))) #this worked... trying to phase out in favor of warning message
      #warning message
      shinyalert::shinyalert(
        title = "Warning: Invalid Input",
        text = paste("The following genes were not found:", paste(missing_genes, collapse = ", ")),
        type = "warning")
    } else {
      # Generate feature plot
      feature_plot <- FeaturePlot(object = seurat_data, features = gene_names,reduction = "umap")
      
      #assign featureplot to global env
      assign("feature_plot", "new", envir = .GlobalEnv) #this doesnt seem to work... try declaring reactive vals at top of server script
      
      #make sure featureplot is not null before rendering plot
      if (is.null(feature_plot)) {
        output$nogenesfound <- renderText("No genes found!")
      } else {
        output$featurePlotUMAP <- renderPlot({ feature_plot })
      }
    }
  })
  
  #pca featureplot logic
  observeEvent(input$plotButton_PCA, {
    req(seuratData(),input$featuresInput) #wait for file input and feature input
    seurat_data <- seuratData()
    if (is.null(seurat_data))
      output$errormessage <- renderText("Null file")
    
    # Split the input string @ commas into individual gene names using regex
    #[[1]] extracts the first element in list of gene names, which gives us a chr vector containing the individaul gene names as separate elements
    gene_names <- strsplit(input$featuresInput, ",\\s*")[[1]] 
    gene_names <- trimws(gene_names)  # Trim leading and trailing spaces from gene names
    
    # Check if any of the requested genes are missing
    missing_genes <- setdiff(gene_names, rownames(seurat_data@assays$RNA@data))
    
    if (length(missing_genes) > 0) {
      # output$nofeaturefound <- renderText(paste("The following genes were not found:", paste(missing_genes, collapse = ", "))) #this worked... trying to phase out in favor of warning message
      #warning message
      shinyalert::shinyalert(
        title = "Warning: Invalid Input",
        text = paste("The following genes were not found:", paste(missing_genes, collapse = ", ")),
        type = "warning")
    } else {
      # Generate feature plot
      feature_plot <- FeaturePlot(object = seurat_data, features = gene_names, reduction = "pca")
      
      #assign featureplot to global env
      # assign("feature_plot", "new", envir = .GlobalEnv) #this doesnt seem to work... try declaring reactive vals at top of server script
      
      #make sure featureplot is not null before rendering plot
      if (is.null(feature_plot)) {
        output$nogenesfound <- renderText("No genes found!")
      } else {
        output$featurePlotPCA <- renderPlot({ feature_plot })
      }
    }
  })
  
  #sampleGene logic
  output$sampleGenes <- renderText({
    seurat_data <- seuratData()
    if (is.null(seurat_data))
      return(NULL)
    # Access the RNA assay
    rna_assay <- seurat_data@assays$RNA
    
    # Get the gene names
    gene_names <- rownames(rna_assay@data)
    
    # Set the random seed for reproducibility (optional)
    set.seed(123)
    
    # Select 10 random genes
    random_genes <- sample(gene_names, 10)
    
    # Return the random genes as a character string
    paste(random_genes, collapse = ", ")
  })
  
  
  ###update dimplot/violinplot dropdown###
  observeEvent(input$seuratFile, {
    req(seuratData()) # wait for file input
    seurat_data <- seuratData()
    
    #fetch metadata colnames
    obj_meta <- seurat_data@meta.data
    
    # Specify the metadata columns you want to check and convert
    columns_to_check <- c("sample", "cytotoxic", "exact_subclonotype_id")
    
    # Find and convert columns that match the specified names
    for (col in columns_to_check) {
      if (col %in% colnames(obj_meta) && !is.factor(obj_meta[[col]])) {
        obj_meta[[col]] <- as.factor(obj_meta[[col]])
      }
    }
    
    #assign column factors to categorical_cols var
    # categorical_cols <- names(obj_meta)[sapply(obj_meta, is.factor)]
    categorical_cols <- names(obj_meta)[sapply(obj_meta, function(col) is.factor(col) && nlevels(col) < 5)]
    
    #reassign dropdown options
    updateSelectInput(session, "variableInput",label = "Select a Field to Split Dim Plot By", choices=categorical_cols)
    
    updateSelectInput(session, "violinInput",label = "Select a Field to Split Violin Plot By", choices = categorical_cols)
    
  })
  
  #new dimplot logic
  observeEvent(input$plotButton_Dim, {
    req(seuratData()) # wait for file input
    seurat_data <- seuratData()
    
    if (is.null(seurat_data)) {
      output$errormessage <- renderText("Null file")
    } else {
      if (!input$splitToggle) {
        # Unsplit dim plot colored by group
        dim_plot <- DimPlot(object = seurat_data, reduction = "umap", group.by = input$variableInput)
        
        if (is.null(dim_plot)) {
          output$nodimplot <- renderText("Bad argument")
        } else {
          output$featurePlotDim <- renderPlot({ dim_plot })
          output$splitToggle <- renderUI({ NULL }) # Hide the splitToggle checkbox
        }
      } else {
        # Split dim plot colored by group
        dim_plot <- DimPlot(object = seurat_data, reduction = "umap", split.by = input$variableInput, group.by = input$variableInput)
        
        if (is.null(dim_plot)) {
          output$nodimplot <- renderText("Bad argument")
        } else {
          output$featurePlotDim <- renderPlot({ dim_plot })
          output$splitToggle <- renderUI({ checkboxInput("splitToggle", "Split Dim Plot", value = TRUE) }) # Show the splitToggle checkbox
        }
      }
    }
  })
  
  ###
  
  #update violinplot dropdown#
  
  
  ### violinplot logic ###
  observeEvent(input$plotButton_violin, {
    req(seuratData(),input$featuresInput_vln,input$violinInput) #wait for all inputs
    seurat_data <- seuratData()
    
    # Split the input string @ commas into individual gene names using regex
    gene_names <- strsplit(input$featuresInput_vln, ",\\s*")[[1]] 
    gene_names <- trimws(gene_names)  # Trim leading and trailing spaces from gene names
    
    # Check if any of the requested genes are missing
    missing_genes_vln <- setdiff(gene_names, rownames(seurat_data@assays$RNA@data))
    
    if (length(missing_genes_vln) > 0) {
      #warning message
      shinyalert::shinyalert(
        title = "Warning: Invalid Input",
        text = paste("The following genes were not found:", paste(missing_genes_vln, collapse = ", ")),
        type = "warning")
    } else {
      # Generate feature plot
      violin_plot <- VlnPlot(object = seurat_data, features = gene_names, split.by =input$violinInput)
      
      #make sure featureplot is not null before rendering plot
      if (is.null(violin_plot)) {
        output$badviolinplot <- renderText("Bad argument")
      } else {
        output$violinplot <- renderPlot({ violin_plot })
      }
    }
  })
} #end of server

#########           SERVER END          ##########


#launch shinyapp
shinyApp(ui = ui, 
         server = server
)