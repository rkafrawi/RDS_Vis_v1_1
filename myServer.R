# load vis source files:

#define server
server <- function(input, output, session) {
  
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
  
  # Render counts table (head)
  output$table <- DT::renderDataTable({
    #retrieve seurat object from reactive expression in above code block
    seurat_data <- seuratData() 
    if (is.null(seurat_data))
      return(NULL)
    # Metadata
    seurat_df <- as.data.frame(seurat_data@meta.data)
    datatable(seurat_df)
  })
  
  ### featureplot logic ###
  #umap featureplot logic
  observeEvent(input$plotButton_UMAP, {
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
      assign("feature_plot", "new", envir = .GlobalEnv) #this doesnt seem to work... try declaring reactive vals at top of server script
      
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
  
  #sample metadata columns logic
  # output$sampleCols <- renderText({
  #   seurat_data <- seuratData()
  #   if (is.null(seurat_data))
  #     return(NULL)
  #   # Access the cols
  #   seurat_cols <- seurat_data@meta.data
  #   
  #   # Get the gene names
  #   seurat_colnames <- names(seurat_cols)
  #   
  #   
  #   # Return the random genes as a character string
  #   paste(seurat_colnames, collapse = ", ")
  # })
  
  
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
    

    #reassign dropdown options (OLD)
    updateSelectInput(session, "variableInput",label = "Select a Field to Split Dimplot By", choices=categorical_cols)

    updateSelectInput(session, "violinInput",label = "Select a Field to Split Dimplot By", choices = categorical_cols)


    
    # Check and convert metadata columns that are not factors ####old way: ERRORS IF THE COLUMNS DONT EXIST#####
    # for (col in columns_to_check) {
    #   if (!is.factor( obj_meta[[col]])) {
    #     obj_meta[[col]] <- as.factor( obj_meta[[col]])
    #   }
    # }
    #categorical colnames only
    # categorical_cols <- names(obj_meta)[sapply(obj_meta, is.factor)]
    #render cat cols
    # output$catcols <- renderText({categorical_cols})
  })
  
  ###dimplot logic### 
  # observeEvent(input$plotButton_Dim, {
  #   req(seuratData()) # wait for file input
  #   seurat_data <- seuratData()
  #   if (is.null(seurat_data)) {
  #     output$errormessage <- renderText("Null file")
  #   } else {
  #     # variable <- switch(input$variableInput,
  #     #                    "Clusters" = "seurat_clusters",
  #     #                    "Treatments" = "sample",
  #     #                    "Clonotype" = "exact_subclonotype_id",
  #     #                    "Cytotoxic" = "cytotoxic")
  # 
  # 
  #     # Generate dim plot
  #     dim_plot <- DimPlot(object = seurat_data, reduction = "umap", split.by=input$variableInput, group.by = input$variableInput) #input$variableInput (use this when autopopulate dropdown)
  #     
  #     # make sure dimplot is not null before rendering plot
  #     if (is.null(dim_plot)) {
  #       output$nodimplot <- renderText("Bad argument")
  #     } else {
  #       output$featurePlotDim <- renderPlot({ dim_plot })
  #     }
  #   }
  # })
  # 
  ##### WIP #####
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
      # violinSelect <- switch(input$violinInput,
      #                    "None" = NULL,
      #                    "Clusters" = "seurat_clusters",
      #                    "Treatments" = "sample")
      
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
  
  ##### WIP #####
  
  
  
  
  ###Volcano plot logic### (7/5/23) incomplete... need sung feedback on what user inputs should be
  # observeEvent(input$plotButton_Volc, {
  #   req(seuratData()) # wait for file input
  #   seurat_data <- seuratData()
  #   if (is.null(seurat_data)) {
  #     output$errormessage <- renderText("Null file")
  #   } else {
  #     # variable <- switch(input$variableInput,
  #     #                    "Clusters" = "seurat_clusters",
  #     #                    "Treatments" = "sample",
  #     #                    "Clonotype" = "exact_subclonotype_id",
  #     #                    "Cytotoxic" = "cytotoxic")
  #     
  #     # Generate volc plot
  #     Volc_plot <- DimPlot(object = seurat_data, reduction = "umap", split.by = input$variableInput)
  #     
  #     # make sure volc is not null before rendering plot
  #     if (is.null( Volc_plot)) {
  #       output$noVolc <- renderText("Bad argument")
  #     } else {
  #       output$Volcplot <- renderPlot({  Volc_plot })
  #     }
  #   }
  # })
  
  # #samplesplit logic
  # output$samplesplit <- renderText({
  #   seurat_data <- seuratData()
  #   if (is.null(seurat_data))
  #     return(NULL)
  #   # Access the metadata
  #   meta <- seurat_data@meta.data
  #   
  #   meta_cols <-names(meta)
  #   
  #   
  #   # Return the params as a character string
  #   paste(meta_cols, collapse = ", ")
  # })
    
  
  # troubleshooting ui x server comms
  
  # Render the input text
  # output$inputz <- renderText({
  #   # input$featuresInput
  #   gene_names
  # })
  #test to see if webapp communicate ui and server
  # output$outputText <- renderPrint({
  #   input$textIn
  # })
  
  ###
}
    
  
    