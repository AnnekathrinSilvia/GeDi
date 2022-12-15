#' GeDi main function
#'
#' @param genesets a dataframe of genesets containing the columns: genesets (the names / ids of the sets) and genes (the genes included in the genesets)
#' @param ppi a Protein-Protein interaction matrix
#' @param alpha a scaling factor in between 0 and 1
#'
#' @return A Shiny app object is returned
#' @export
#' @import visNetwork
#' @import shiny
#' @import shinyBS
#' @import fontawesome
#' @importFrom rintrojs introjs
#' @importFrom utils read.delim
#' @importFrom bs4Dash bs4DashPage bs4DashNavbar box bs4DashBrand bs4DashBody bs4Card bs4DashSidebar bs4SidebarMenu bs4SidebarMenuItem bs4TabItem bs4TabItems tabBox
#' @importFrom shinycssloaders withSpinner
#'
#' @examples
#' \dontrun{
#' GeDi()
#' }
GeDi <- function(genesets = NULL,
                 ppi = NULL,
                 alpha = 1) {
  options(spinner.type = 6, spinner.color = "#0092AC")
  #on.exit(options(oopt))

  usage_mode <- "shiny_mode"

  # UI definition -----------------------------------------------------------

  # dashpage definition -----------------------------------------------------
  gedi_ui <- bs4DashPage(
    title = "GeDi",
    dark = NULL,
    # navbar definition -------------------------------------------------------
    header = bs4DashNavbar(
      tagList(tags$code(tags$h3("GeDi"))),
      tags$span(style = "display:inline-block; width: 30%"),
      tagList(
        shinyWidgets::dropdownButton(
          inputId = "ddbtn_docs",
          circle = FALSE,
          status = "info",
          icon = icon("book"),
          width = "300px",
          size = "xs",
          right = TRUE,
          tooltip = shinyWidgets::tooltipOptions(title = "More documentation on GeDi"),
          tags$h5("Documentation"),
          actionButton(
            inputId = "btn_first_help",
            icon = icon("circle-question"),
            label = "First Help",
            style = .actionButtonStyle
          )
        ),
        shinyWidgets::dropdownButton(
          inputId = "ddbtn_info",
          circle = FALSE,
          status = "info",
          icon = icon("info"),
          width = "300px",
          size = "xs",
          right = TRUE,
          tooltip = shinyWidgets::tooltipOptions(title = "More info on GeDi and on the current session"),
          tags$h5("Additional information"),
          actionButton(
            inputId = "btn_info_session",
            icon = icon("circle-info"),
            label = "About this session",
            style = .actionButtonStyle
          ),
          actionButton(
            inputId = "btn_info_genetonic",
            icon = icon("heart"),
            label = "About GeDi",
            style = .actionButtonStyle
          )
        )
      )
    ),
    # ),
    # title = bs4DashBrand(
    #   title = HTML("<small>GeDi</small>"),
    #   #href = "https://bioconductor.org/packages/GeneTonic",
    #   # color = "info",
    #   #image = "GeneTonic/GeneTonic.png"),
    #   skin = "dark",
    #   status = "gray-dark",
    #   border = FALSE,
    #   controlbarIcon = icon("cogs"),
    #   fixed = TRUE
    # ),

    # sidebar definition ------------------------------------------------------
    sidebar = bs4DashSidebar(
      id = "sidebar",
      title = HTML("<small>GeDi</small>"),
      #src = "GeneTonic/GeneTonic.png",
      skin = "dark",
      status = "info",
      brandColor = NULL,
      #url = "https://bioconductor.org/packages/GeneTonic",
      collapsed = TRUE,
      elevation = 1,
      opacity = 0.8,
      bs4SidebarMenu(
        bs4SidebarMenuItem("Welcome!",
                           tabName = "tab_welcome",
                           icon = icon("house")
                           ),
        bs4SidebarMenuItem(
          "Data Upload",
          tabName = "tab_data_upload",
          icon = icon("square-share-nodes")
          ),
        bs4SidebarMenuItem(
          "Distance Scores",
          tabName = "tab_scores",
          icon = icon("diagram-project")
          ),
        bs4SidebarMenuItem(
          "Graph",
          tabName = "tab_graph",
          icon = icon("square-share-nodes")
          )
      )
    ),

    # body definition ---------------------------------------------------------
    body = bs4DashBody(
      rintrojs::introjsUI(),
      tags$head(tags$style(
        HTML(
          ".shiny-output-error-validation {
            font-size: 15px;
            color: forestgreen;
            text-align: center;
            }
            #myScrollBox{
              overflow-y: scroll;
              .dataTables_wrapper{
                overflow-x: scroll;
              }
            }
          .tooltip > .tooltip-inner {
                 width: 1000px;
                 color: white;
                 background-color: black;
                 }"
          )
        )
      ),
      tags$head(
        tags$style(
          ".biocdlbutton{background-color:#0092AC;} .biocdlbutton{color: #ffffff;}"
        )
      ),
      tags$script(
        HTML(
          "$(function(){
      $(document).keyup(function(e) {
      if (e.which == 17) {
        $('#bookmarker').click()
      }
      });
      })"
        )
      ),
      bs4TabItems(
        # ui panel welcome ---------------------------------------------------
        bs4TabItem(tabName = "tab_welcome",
                   tagList(
                     fluidRow(column(width = 11),
                              column(
                                width = 1,
                                actionButton(
                                  "tour_welcome",
                                  label = "",
                                  icon = icon("circle-question"),
                                  style = .tourButtonStyle
                                  )
                                )
                              ),
                     uiOutput("ui_panel_welcome")
                     )
                   ),
        # ui panel data upload -----------------------------------------------
        bs4TabItem(
          tabName = "tab_data_upload",
          tagList(
            fluidRow(
              column(width = 11),
              column(
                width = 1,
                br(),
                actionButton(
                  "tour_data_upload",
                  label = "",
                  icon = icon("circle-question"),
                  style = .tourButtonStyle
                ),
                shinyBS::bsTooltip(
                  id = "tour_data_upload",
                  title = "Click me to start a tour of this section!",
                  placement = "top",
                  trigger = "hover",
                  options = list(container = "body")
                ),
                style = "float:right"
              )
            ),
            uiOutput("ui_panel_data_upload"),
            uiOutput("ui_panel_specify_species"),
            uiOutput("ui_panel_download_ppi")
          )
        ),
        # ui panel scores -------------------------------------------------
        bs4TabItem(tabName = "tab_scores",
                   tagList(
                     fluidRow(
                       column(width = 11),
                       column(
                         width = 1,
                         br(),
                         actionButton(
                           "tour_scoring",
                           label = "",
                           icon = icon("circle-question"),
                           style = .tourButtonStyle
                         ),
                         shinyBS::bsTooltip(
                           id = "tour_scoring",
                           title = "Click me to start a tour of this section!",
                           placement = "top",
                           trigger = "hover",
                           options = list(container = "body")
                         ),
                         style = "float:right"
                       )
                     ),
                     uiOutput("ui_panel_scores")
                     )
                   ),
        # ui panel graph ---------------------------------------------------
        bs4TabItem(tabName = "tab_graph",
                   tagList(
                     fluidRow(
                       column(width = 11),
                       column(
                         width = 1,
                         br(),
                         actionButton(
                           "tour_graph",
                           label = "",
                           icon = icon("circle-question"),
                           style = .tourButtonStyle
                         ),
                         shinyBS::bsTooltip(
                           id = "tour_graph",
                           title = "Click me to start a tour of this section!",
                           placement = "top",
                           trigger = "hover",
                           options = list(container = "body")
                         ),
                         style = "float:right"
                       )
                     ),
                     uiOutput("ui_panel_graph")
                     )
                   )
      )
    )
    # ),
    # # controlbar definition ------------------------------------------------
    # controlbar = bs4Dash::bs4DashControlbar(
    #   collapsed = TRUE,
    #   uiOutput("ui_controlbar")
    # ),
    #
    # # footer definition -------------------------------------------------------
    # footer = bs4DashFooter(
    #   left = GeneTonic_footer,
    #   right = NULL
    # )
  )
  # server deifnition ---------------------------------------------------------
  gedi_server <- function(input, output, session) {
    # initializing reactives --------------------------------------------------
    reactive_values <- reactiveValues()
    reactive_values$genesets <- NULL
    reactive_values$gs_names <- NULL
    reactive_values$genes <- NULL
    reactive_values$species <- NULL
    reactive_values$ppi <- NULL
    reactive_values$scores <- NULL
    reactive_values$seeds <- NULL
    reactive_values$cluster <- NULL


    # panel Welcome ----------------------------------------------------------
    output$ui_welcome <- renderUI({
      tagList(fluidRow(column(width = 12)))
    })


    # panel Data Upload ------------------------------------------------------
    output$ui_panel_data_upload <- renderUI({
      tagList(
        box(
          width = 12,
          title = "Step 1",
          status = "danger",
          solidHeader = TRUE,
          h2("Upload your Genesets input data"),
          fluidRow(
            column(
              width = 6,
              uiOutput("upload_genesets")
              # br(),
              # "... or you can also ",
              # actionButton("btn_loaddemo", "Load the demo airway data",
              #              icon = icon("play-circle"),
              #              class = "btn btn-info"
              # ), br(), p()),
            ),
            column(
              width = 6,
              box(
                id = "Genests_preview",
                width = NULL,
                title = "Genesets preview",
                status = "info",
                solidHeader = TRUE,
                collapsible = TRUE,
                collapsed = TRUE,
                fluidRow(column(
                  width = 12,
                  offset = 0.5,
                  DT::dataTableOutput("dt_genesets")
                  )
                )
              )
            )
          )
        )
      )
    })

    output$upload_genesets <- renderUI({
      return(
        fileInput(
          inputId = "uploadgenesetfile",
          label = "Upload a geneset file",
          accept = c(
            "text/csv",
            "text/comma-separated-values",
            "text/tab-separated-values",
            "text/plain",
            ".csv",
            ".tsv"
          ),
          multiple = FALSE
        )
      )
    })

    readGenesets <- reactive({
      if (is.null(input$uploadgenesetfile)) {
        return(NULL)
      }

      guessed_sep <-
        sepguesser(input$uploadgenesetfile$datapath)
      genesets <-
        utils::read.delim(
          input$uploadgenesetfile$datapath,
          header = TRUE,
          as.is = TRUE,
          sep = guessed_sep,
          quote = "",
          check.names = FALSE
        )
      return(genesets)
    })

    output$dt_genesets <- DT::renderDataTable({
      validate(need(!(
        is.null(reactive_values$genesets)
        ),
      message = "Please upload a text file via the button on the left.")
      )

      DT::datatable(reactive_values$genesets,
                    options = list(scrollX = TRUE, scrollY = "400px"))
    })


    output$ui_panel_specify_species <- renderUI({
      if (is.null(reactive_values$genesets)) {
        return(NULL)
      }
      box(
        width = 12,
        title = "Step 2",
        status = "warning",
        solidHeader = TRUE,
        tagList(
          h2("Select the species of your data"),
          fluidRow(column(
            width = 6,
            uiOutput("ui_specify_species")
          ),
          column(
            width = 6,
            uiOutput("ui_specify_opt_params_ppi_download")
            )
          )
        )
      )
    })

    output$ui_specify_species <- renderUI({
      if (is.null(reactive_values$genesets)) {
        return(NULL)
      }
      fluidRow(column(
        width = 12,
        textInput(
          "species",
          label = "Please specify the species of your data
                          (e.g. Homo Sapiens, Mus musculus, etc.)",
          width = '400px'
          )
        )
      )
    })

    output$ui_specify_opt_params_ppi_download <- renderUI({
      if (is.null(reactive_values$genesets)) {
        return(NULL)
      }
      box(
        id = "Optional_parameters",
        width = 12,
        title = "Optional Parameters",
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = TRUE,
        tagList(fluidRow(column(
          width = 6,
          fluidRow(column(
            width = 12,
            textInput("stringVersion",
                      label = "Specify the StringDb version to use",
                      value = "11.5"),
            textInput(
              "scoreThresholdString",
              label = "Specify the score threshold",
              value = "0.00"
            ),
            textInput(
              "inputDirectoryString",
              label = "Specify the input directory",
              value = ""
              )
            )
            )
          )
          )
        )
      )
    })

    output$ui_panel_download_ppi <- renderUI({
      if (input$species == "" ||
          is.na(input$species) || is.null(input$species)) {
        return(NULL)
      }
      reactive_values$species <- input$species
      box(
        width = 12,
        title = "Step 3",
        status = "success",
        solidHeader = TRUE,
        tagList(
          h2("Download the PPI matrix from STRING"),
          fluidRow(column(
            width = 6,
            actionButton(
              "download_ppi",
              label = "Download PPI matrix",
              icon = icon("download"),
              style = .actionButtonStyle
            )
          ),
          column(
            width = 6,
            box(
              id = "PPI_preview",
              width = NULL,
              title = "PPI preview",
              status = "info",
              solidHeader = TRUE,
              collapsible = TRUE,
              collapsed = TRUE,
              fluidRow(column(
                width = 12,
                offset = 0.5,
                DT::dataTableOutput("dt_ppi")
                )
                )
              )
            )
          )
        )
      )
    })

    output$dt_ppi <- DT::renderDataTable({
      validate(
        need(!(is.null(reactive_values$ppi)),
        message = "Please download a Protein-Protein Interaction (PPI) matrix via the button on the left.")
      )
      DT::datatable(reactive_values$ppi,
                    options = list(scrollX = TRUE, scrollY = "400px"))
    })

    # panel Scores ----------------------------------------------------
    output$ui_panel_scores <- renderUI({
      tagList(
        box(
          id = "distance_calc_box",
          width = 12,
          title = "Calculate distance scores for you Genesets",
          status = "info",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = FALSE,
          h2("Select the distance score"),
          fluidRow(
            column(
              width = 6,
              radioButtons(
                inputId = "scoringmethod",
                label = "Select the scoring method for your data",
                choices = c("Meet-Min", "Kappa", "PMM", "Jaccard"),
                selected = character(0)
              )
            ),
            column(width = 6,
                   uiOutput("ui_score_data"))
          ),
          fluidRow(column(width = 12))
        ),
        fluidRow(column(
          width = 12,
          bs4Dash::bs4Card(
            id = "tabcard_scores",
            title = "Geneset Distance Scores",
            elevation = 1,
            width = 12,
            closable = FALSE,
            bs4Dash::tabsetPanel(
              id = "tabsetpanel_scores",
              type = "tabs",
              selected = "Distance Scores Heatmap",
              side = "right",
              tabPanel(title = "Distance Scores Heatmap",
                       withSpinner(plotOutput(
                         "scores_heatmap"))
                       ),
              tabPanel(title = "Distance Scores Graph",
                       fluidRow(
                         column(
                           width = 2,
                           sliderInput(
                             "similarityScores",
                             "Distance Threshold:",
                             min = 0,
                             max = 1,
                             value = 0.3,
                             step = 0.05
                           )
                         ),
                         column(width = 9,
                                withSpinner(
                                  visNetworkOutput("scores_Network")
                                  )
                                )
                         )
                       )
              )
            )
          )
        )
      )
    })

    # TODO: Workaround with empty string, check to fix that later in another way
    output$ui_score_data <- renderUI({
      if (is.null(input$scoringmethod) || input$scoringmethod == "") {
        return(NULL)
      }
      fluidRow(
        column(
          width = 12,
          strong("Now you can score your data"),
          br(),
          "Attention: If you have many genesets to score,
           this operation may take some time",
          br(),
          actionButton("score_data",
                       label = "Score the Genesets",
                       style = .actionButtonStyle)
        )
      )
    })

    output$scores_heatmap <- renderPlot({
      validate(need(!(is.null(reactive_values$scores)),
      message = "Please score you genesets first in the above box"))
      distance_heatmap(reactive_values$scores,
                       chars_limit = 20)
    })

    reactive_values$scores_graph <- reactive({
      adj <- getAdjacencyMatrix(reactive_values$scores,
                                input$similarityScores)
      g <- buildGraph(adj)
      return(g)
    })

    output$scores_Network <- renderVisNetwork({
      validate(need(!(is.null(reactive_values$scores)),
      message = "Please score you genesets first in the above box"))


      if (!any(igraph::get.edgelist(reactive_values$scores_graph()) != 0)) {
        showNotification(
          "Please select a larger distance threshold as currently no nodes are connected and the graph cannot be properly rendered.",
          type = "warning"
        )
      } else{
        visNetwork::visIgraph(reactive_values$scores_graph()) %>%
          visOptions(
            highlightNearest = list(
              enabled = TRUE,
              degree = 1,
              hover = TRUE
            ),
            nodesIdSelection = TRUE
          ) %>%
          visExport(name = "backbone_network",
                    type = "png",
                    label = "Save backbone graph")
      }
    })


    # panel Graph ------------------------------------------------------

    output$ui_panel_graph <- renderUI({
      validate(need(!is.null(reactive_values$scores),
                    message = "Please score you genesets first in the Scores panel."))

      tagList(
        box(
          id = "clustering_param_box",
          width = 12,
          title = "Set the parameters for Clustering",
          status = "info",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = FALSE,
          fluidRow(
            column(
              width = 6,
              sliderInput(
                inputId = "simThreshold",
                label = "Similarity Threshold",
                min = 0,
                max = 1,
                value = 0.3,
                step = 0.1
              ),
              sliderInput(
                inputId = "memThreshold",
                label = "Membership Threshold",
                min = 0,
                max = 1,
                value = 0.5,
                step = 0.1
              ),
              sliderInput(
                inputId = "clustThreshold",
                label = "Clustering Threshold",
                min = 0,
                max = 1,
                value = 0.5,
                step = 0.1
              )
            ),
            column(width = 6,
                   uiOutput("ui_cluster"))
          ),
          fluidRow(column(width = 12))
        ),
        fluidRow(
          column(
            width = 12,
            bs4Dash::bs4Card(
              id = "tabcard_cluster",
              title = "Geneset Cluster Graph",
              elevation = 1,
              width = 12,
              closable = TRUE,
              bs4Dash::tabsetPanel(
                id = "tabsetpanel_cluster",
                type = "tabs",
                selected = "Geneset Graph",
                side = "right",
                tabPanel(title = "Geneset Graph",
                         withSpinner(
                           visNetworkOutput("cluster_Network",
                                            height = "700px",
                                            width = "100%")
                         )
                ),
                tabPanel(title = "Cluster-Geneset Bipartite Graph",
                         withSpinner(
                           visNetworkOutput("cluster_geneset_bipartite_Network",
                                            height = "700px",
                                            width = "100%")
                           )
                         )
                )
              )
            )
          ),
        fluidRow(
          bs4Dash::bs4Card(
            width = 12,
            id = "card_clustering",
            title = "Clustering graph summaries",
            status = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            closable = FALSE,
            fluidRow(column(
              width = 12,
              DT::dataTableOutput("dt_cluster")
              )
              )
          )
        )
      )
    })

    output$ui_cluster <- renderUI({
      fluidRow(
        column(
          width = 12,
          strong("Now you can cluster your data"),
          br(),
          "Attention: If you have many Genesets to cluster,
           this operation may take some time",
          br(),
          actionButton("cluster_data",
                       label = "Cluster the Genesets",
                       style = .actionButtonStyle)
        )
      )
    })

    output$cluster_Network <- renderVisNetwork({
      validate(need(!(is.null(reactive_values$cluster)),
      message = "Please cluster you genesets first in the above box"))


      if (!any(igraph::get.edgelist(reactive_values$cluster_graph()) != 0)) {
        showNotification(
          "It seems like you don't have any clusters. Please adapt the similarity Threshold above and re-run the Clustering.",
          type = "warning"
        )
      } else{
        visNetwork::visIgraph(reactive_values$cluster_graph()) %>%
          visOptions(
            highlightNearest = list(
              enabled = TRUE,
              degree = 1,
              hover = TRUE
            ),
            nodesIdSelection = TRUE
          ) %>%
          visExport(name = "cluster_network",
                    type = "png",
                    label = "Save Cluster graph")
      }
    })

    reactive_values$cluster_graph <- reactive({
      adj <- getClusterAdjacencyMatrix(reactive_values$cluster,
                                       reactive_values$gs_names)
      g <- buildGraph(adj)
      return(g)
    })



    output$cluster_geneset_bipartite_Network <- renderVisNetwork({
      validate(need(!(is.null(reactive_values$cluster)),
                    message = "Please cluster you genesets first in the above box"))

        visNetwork::visIgraph(reactive_values$bipartite_graph()) %>%
          visIgraphLayout(layout = "layout_as_bipartite") %>%
          visOptions(
            highlightNearest = list(
              enabled = TRUE,
              degree = 1,
              hover = TRUE
            ),
            nodesIdSelection = TRUE
          ) %>%
          visExport(name = "bipartite_network",
                    type = "png",
                    label = "Save Cluster-Geneset bipartite graph")
    })

    reactive_values$bipartite_graph <- reactive({
      g <- getBipartiteGraph(reactive_values$cluster,
                             reactive_values$gs_names,
                             reactive_values$genes
                             )
      #g <- add_layout_(as_bipartite(g))
      return(g)
    })


    output$dt_cluster <- DT::renderDataTable({
      validate(need(!(is.null(reactive_values$cluster)),
        message = "It seems like the data has not yet been clustered. Please cluster your data first with the box above.")
      )
      dt_cluster <- getClusterDatatable(reactive_values$cluster,
                                        reactive_values$gs_names)
      DT::datatable(dt_cluster,
                    options = list(scrollX = TRUE, scrollY = "400px"))
    })

    # controlbar --------------------------------------------------------------

    # output$ui_controlbar <- renderUI({
    #   validate(
    #     need(!is.null(reactive_values$gtl) ,
    #          message = "\n\n\nPlease provide a GeneTonicList or its components to display the content of the control bar.")
    #   )
    #
    #   # message(nrow(reactive_values$res_enrich))
    #
    #   tagList(
    #     numericInput(
    #       inputId = "de_fdr",
    #       label = "False Discovery Rate (FDR) for DE",
    #       value = 0.05, min = 0.0001, max = 1, step = 0.01
    #     ),
    #     numericInput(
    #       inputId = "n_genesets",
    #       label = "Number of genesets",
    #       value = 15, min = 1, max = nrow(reactive_values$res_enrich)
    #     ),
    #     selectInput("exp_condition",
    #                 label = "Group/color by: ",
    #                 choices = c(NULL, names(colData(reactive_values$dds))), selected = NULL, multiple = TRUE
    #     ),
    #     colourInput("col", "Select colour for volcano plot", "#1a81c2",
    #                 returnName = TRUE,
    #                 allowTransparent = TRUE
    #     ),
    #     checkboxInput("labels", label = "Display all labels", value = FALSE)
    #   )
    # })

    # outputOptions(output, "ui_controlbar", suspendWhenHidden = FALSE)

    # Observers ------------------------------------------------------------------------------------
    observeEvent(input$uploadgenesetfile, {
      reactive_values$genesets <- readGenesets()
      reactive_values$gs_names <- reactive_values$genesets$Geneset
      tryCatch(
        expr = {
          reactive_values$genes <- getGenes(reactive_values$genesets)
        },
        error = function(cond) {
          showNotification(
            "It seems like your data does not have a column named 'Genes'. Please check your data and try to upload it again.",
            type = "error"
          )
          reactive_values$genesets <- NULL
          reactive_values$gs_names <- NULL
          reactive_values$genes <- NULL
        }
      )
    })

    observeEvent(input$download_ppi, {
      validate(need(!(input$species == ""),
                    message = "Please specify the species of your data"))
      reactive_values$species <- input$species
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "Downloading PPI from STRINGdb", value = 0)
      progress$inc(1 / 12, detail = "Get species ID")
      id <- getId(reactive_values$species)
      progress$inc(1 / 12, detail = "Validate species ID")
      validate(need(!is.na(id),
          message = "We could not find your specified species.
                   Please check the spelling and try again."
        )
      )
      progress$inc(1 / 12, detail = "Get species specific STRINGdb")
      stringdb <- getStringDB(as.numeric(id))
      stringdb

      progress$inc(4 / 12, detail = "Get Annotation information")
      anno_df <- getAnnotation(stringdb)

      progress$inc(1 / 12, detail = "Download PPI")
      reactive_values$ppi <- getPPI(reactive_values$genes,
                                    string_db = stringdb,
                                    anno_df = anno_df)
      progress$inc(4 / 12, detail = "Done")
    })

    observeEvent(input$score_data, {
      if (input$scoringmethod == "" ||
          (is.null(input$scoringmethod))) {
        showNotification(
          "It seems like you did not select a scoring method. Please select a scoring method on the left.",
          type = "error"
        )
      }

      if ((length(reactive_values$genes) == 0  ||
           (is.null(reactive_values$genes)))) {
        showNotification(
          "It seems like the file you've provided does not
             contain any genesets. Please check you input and retry.",
          type = "error"
        )
      }

      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "Scoring your genesets", value = 0)

      if (input$scoringmethod == "Meet-Min") {
        scores <- getMeetMinMatrix(reactive_values$genes, progress)
      } else if (input$scoringmethod == "Kappa") {
        scores <- getKappaMatrix(reactive_values$genes, progress)
      } else if (input$scoringmethod == "PMM") {
        # TODO: Handle alpha as possible input
        if (is.null(reactive_values$ppi)) {
          showNotification(
            "It seems like you have not downloaded a PPI matrix. Please return to the Data Upload panel and download the respective PPI.",
            type = "error"
          )
          scores <- NULL
        } else{
          scores <- getpMMMatrix(reactive_values$genes,
                                 reactive_values$ppi,
                                 progress = progress)
        }
      } else if (input$scoringmethod == "Jaccard") {
        scores <- getJaccardMatrix(reactive_values$genes,
                                   progress = progress)
      }

      if (is.null(scores)) {
        showNotification(
          "It seems like something went wrong while scoring your data. Most likely this is due to the genesets not containg genes.",
          type = "error"
        )
      } else{
        rownames(scores) <- colnames(scores) <- reactive_values$gs_names
        reactive_values$scores <- scores
        bs4Dash::updateBox("distance_calc_box", action = "toggle")
      }
    })

    observeEvent(input$cluster_data, {
      if (is.null(reactive_values$scores)) {
        showNotification(
          "It seems like you did not score your Genesets yet. Please go back to
          the Scores panel and select a score of your choice.",
          type = "error"
        )
      }
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "Finding initial seeds", value = 0)

      seeds <- seedFinding(reactive_values$scores,
                           input$simThreshold,
                           input$memThreshold)

      progress$inc(0.5, detail = "Found initial seeds. Now clustering the data.")

      cluster <- clustering(seeds, input$clustThreshold)

      progress$inc(0.5, detail = "Finished clustering the data.")

      if (is.null(seeds)) {
        showNotification(
          "It seems like something went wrong while clustering your data. Most likely this is due to the genesets not containg genes.",
          type = "error"
        )
      } else{
        reactive_values$seeds <- seeds
        reactive_values$cluster <- cluster
        bs4Dash::updateBox("clustering_param_box", action = "toggle")
      }
    })


    # Tour Observers -------------------------------------------------------------------------
    observeEvent(input$tour_welcome, {
      tour <- utils::read.delim(
        system.file("extdata",
                    "intro_welcome.txt",
                    package = "GeDi"),
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_data_upload, {
      tour <- utils::read.delim(
        system.file("extdata",
                    "intro_data_upload.txt",
                    package = "GeDi"),
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_scoring, {
      tour <- utils::read.delim(
        system.file("extdata",
                    "intro_scoring.txt",
                    package = "GeDi"),
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_graph, {
      tour <- utils::read.delim(
        system.file("extdata",
                    "intro_graph.txt",
                    package = "GeDi"),
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session, options = list(steps = tour))
    })

  }
  shinyApp(ui = gedi_ui, server = gedi_server)
}
