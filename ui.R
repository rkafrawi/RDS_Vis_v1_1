
#########           UI START          ##########

# fully defined ui
dashboardPage(
  
  # webapp layout #
  dashboardHeader(title = "RDS Vis v1"),
  menu_bar <- dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Data Summary", tabName="input_summ", icon=icon("table")),
      menuItem("Feature Plot", tabName = "input", icon = icon("chart-bar")),
      menuItem("Dim Plot", tabName = "dim", icon = icon("cubes")),
      menuItem("Violin Plot", tabName = "violin", icon = icon("record-vinyl")),
      # menuItem("Volcano Plot", tabName = "volcano", icon = icon("fire")),
      menuItem("Help Page", tabName = "help", icon = icon("question"))
    )
  ),
  dashboardBody(
    tabItems(
      tab_home <- tabItem(tabName = "home",
                          h2("Home Page"),
                          br(),
                          br(),
                          p("Welcome to the RDS Vis Tool!"),
                          br(),
                          p("This web app is designed for the visualization and exploration of single-cell RNA-seq data contained in Seurat objects. 
                             It will provide various interactive plots and features to help you analyze and gain insights from your data. 
                             A cloud build of this app is being developed on DNAnexus, which will circumvent the instance size restrictions of shinyapps.io's base plan."),
                          br(),
                          p("This web app currently consists of four pages of note: a data summary page, a feature plot page, a dim plot page, and a violin plot page. 
                            As a general rule of thumb, these pages will not visualize any data if no .rds file has been provided by the user in the data summary page.
                            Note also that the drop downs on the dim plot and violin plot pages will not display updated categories until a file has been provided by the user."),
                          br(),
                          #embed href in paragraph sentence
                          p("Click ", a("here", target="_blank", href="https://raw.githubusercontent.com/rkafrawi/RDS_Vis_v1_1/main/RDS_Vis_1.1.pdf", .noWS = "outside"), " for an in depth user guide.", .noWS = c("after-begin", "before-end")),
                          
                          # Add any additional content or UI elements here
      ),
      tab_input <- tabItem(tabName = "input_summ",
                           h2("Data Summary Page"),
                           br(),
                           #file upload
                           p("Click on the browse button and select the .rds file you wish to generate a feature plot for (Max 10 GB). 
                                     Once your upload has completed, some sample genes from the Seurat object will populate below."),
                           fileInput("seuratFile", 
                                     "Upload RDS File",
                                     accept = ".rds"),
                           br(),
                           p("Table of Seurat object metadata populates below:"),
                           dataTableOutput("table")
      ),
      tab_inputfeatures <- tabItem(tabName = "input",
                                   h2("Feature Plot Page"),
                                   br(),
                                   
                                   br(),
                                   p("Try some of these sample genes:"),
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
                         # selectInput("variableInput", label = NULL,
                         #             choices = c("Clusters", "Treatments", "Clonotype", "Cytotoxic" )),
                         selectInput("variableInput", label = "Select a Group:",
                                     choices = "Input File For Dropdown Options"),
                         
                         br(),
                         
                         br(),
                         actionButton("plotButton_Dim","Generate DimPlot"),
                         br(),
                         plotOutput("featurePlotDim", height = 600, width = 900)
                         
      ),
      tab_violin <- tabItem(tabName = "violin",
                            h2("Violin Plot Page"),
                            # Add violin plot page content here
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
                              title = "Is there any way to view tabular representations of the Seurat object?",
                              "Currently, the server does not have prebuilt logic to generate tables of the Seurat object's metadata outside of the summary table on the Data Summary page."
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
          src ="https://raw.githubusercontent.com/rkafrawi/RDS_Vis_v1_1/main/bxMD_logo.jpeg",
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
