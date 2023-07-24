# load tab source files:
# local 
# list.files(path = "C:/Users/rkafr/Desktop/dnanexus_app_v1/tabs") %>%
#   str_subset("\\.R") %>%
#   map(~source(paste0("C:/Users/rkafr/Desktop/dnanexus_app_v1/tabs/", .x)))

#web server
# source_files <- list.files(path = "./dnanexus_app_v1/tabs", pattern = "\\.R$", full.names = TRUE)
# 
# lapply(source_files, source)


# compressed ui
# ui <- dashboardPage(
#   
#   # webapp layout #
#   dashboardHeader(title = "RDS Vis v1"),
#   menu_bar,#from list.files() fxn
#   dashboardBody(
#     tabItems(
#       tab_home,
#       tab_inputfeatures,
#       # tab_input,
#       # tab_feature,
#       # tab_volcano,
#       # tab_violin,
#       # tab_dim,
#       tab_help
#     ) #everything within tabItems() was fetched from the list.files() fxn
#   )
# )

# fully defined ui
ui <- dashboardPage(
  
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
                             Currently, this webapp does not support use of more than 1 GB of memory due to instance size restrictions of shinyapps.io's base plan. 
                             As such, while the local version that we have been testing does not require the user to refresh between .rds file uploads, the server version you are using in some cases will do so.
                             This webapp in its current state only exists to give users a sense of the page layout and to give developers insight into user requests.
                             This memory restriction will no longer exist once the webapp is deployed on DNAnexus, which will lead to a more smooth and streamlined user experience."),
                          br(),
                          
                          #embed href in paragraph sentence
                          p("Click ", a("here", target="_blank", href="RDS_Vis_1.1.pdf", .noWS = "outside"), " for an in depth user guide.", .noWS = c("after-begin", "before-end")),
                          
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
                                   
                                   # 
                                   # tabPanel("UMAP", plotOutput("featurePlot")),
                                   # br(),
                                   # textOutput("nofeaturefound"),
                                   # br()
                                   
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
                         
                         
                         # i think im going to go for the tab approach. put the dropdown within the subtab.
                         # tabsetPanel(
                         #   tabPanel("Unsplit Dimplot",
                         #            br(),
                         #            selectInput("variableInput_US", label = "Select",
                         #                        choices = "Input File For Dropdown Options"),
                         #            br(),
                         #            actionButton("plotButton_Dim_US","Generate DimPlot"),
                         #            br(),
                         #            plotOutput("featurePlotDim_US", height = 600, width = 900)
                         # 
                         # 
                         #   ),
                         #   tabPanel("Split Dimplot",
                         #            p("Select a Group to split the Dimplot By:"),
                         #            selectInput("variableInput_S", label = "Select",
                         #                        choices = "Input File For Dropdown Options"),
                         # 
                         #            br(),
                         #            actionButton("plotButton_Dim_S","Generate DimPlot"),
                         #            br(),
                         #            plotOutput("featurePlotDim_S", height = 600, width = 900)
                         #   
                         # )),
                         br(),
                         # selectInput("variableInput", label = NULL,
                         #             choices = c("Clusters", "Treatments", "Clonotype", "Cytotoxic" )),
                         selectInput("variableInput", label = "Select a Group:",
                                     choices = "Input File For Dropdown Options"),
                   
                         br(),
                         # p("Here are the available categorical columns that the Dimplot can be split by:"),
                         # br(),
                         # div(verbatimTextOutput("catcols"), class = "full-width-cols"),


                         # p("Metadata split.by params:"),
                         # br(),
                         # verbatimTextOutput("samplesplit"),
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
      
                         # Add any additional content or UI elements here
      # tab_volcano <- tabItem(tabName = "volcano",
      #                        h2("Volcano Plot Page"),
      #                        br(),
      #                        p("Try some of these sample genes:"),
      #                        br(),
      #                        verbatimTextOutput("sampleGenes"),
      #                        br(),
      #                        p("Enter gene name:"),
      #                        textInput("geneVolc", "Gene Name"),
      #                        br(),
      #                        p("Select a variable:"),
      #                        selectInput("volcanoInput", label = NULL,
      #                                    choices = c("Clusters", "Treatments", "Responder States")),
      #                        br(),
      #                        actionButton("plotButton_Volc","Generate DimPlot"),
      #                        br(),
      #                        plotOutput("Volcplot", height = 600, width = 900)
      #                        # Add volcano plot page content here
      # ),
      
      tab_help <- tabItem(tabName = "help",
                          h2("Help Page"),
                          br(),
                          p("Below are a list of potential FAQs (Subject to change)."),
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
                              "Currently, the server does have prebuilt logic to generate tables of the Seurat object's metadata. It is simply not present in the current version of the user interface to prevent visually bogging down the webapp and general ease of use. "
                            )
                          ),
                          br(),
                          br(),
                          br(),
                          #placeholder to make footer show on help page
                          div(id = "placeholder",p("~ More FAQs in Progress! ~"))
      )
    ), 
    div(class = "hr-container",
        hr(class = "hr-line")),
    # Footer
    tags$footer(
      class = "footer",
      p("Developed for the Belfer Center Research Team by the BxMD Group at DFCI.")
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